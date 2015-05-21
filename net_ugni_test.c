/*
  This file is part of SystemConfidence.

  Copyright (C) 2012, UT-Battelle, LLC.

  This product includes software produced by UT-Battelle, LLC under Contract No. 
  DE-AC05-00OR22725 with the Department of Energy. 

  This program is free software; you can redistribute it and/or modify
  it under the terms of the New BSD 3-clause software license (LICENSE). 
  
  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the 
  LICENSE for more details.

  For more information please contact the SystemConfidence developers at: 
  systemconfidence-info@googlegroups.com

*/

/**
 * \brief Tests a single (simultaneous, bidirectional) exchange, but extracts
 * both one-sided and pairwise variability.
 *
 * Pros: 
 * - Assesses both one-sided and pairwise variability with minimal averaging
 * - Provides a Least Upper Bound of the network's minimum latency
 * - Quantifies network topology effects 
 * - Provides a baseline minimum for comparison
 *
 * Cons:
 * - Requires additional storage
 */

#include <unistd.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <stdio.h>
#include <errno.h>

#include "config.h"
#include "orbtimer.h"
#include "comm.h"
#include "tests.h"
#include "measurement.h"
#include "gni_pub.h"

#ifdef SHMEM
	#include <mpp/shmem.h>
#else
	#include <mpi.h>
#endif

/* number of network latency histograms */
#define NET_LEN 9

/* histograms for network latency measurements */
enum net_vars {
	/* timer overhead */
	timer,
	/* local communication */
	onNodeOnesided, onNodePairwise,
	onNodeOnesidedMinimum, onNodePairwiseMinimum,
	/* remote communication */
	offNodeOnesided, offNodePairwise,
	offNodeOnesidedMinimum, offNodePairwiseMinimum
};

char *net_labels[] = {
	/* timer overhead */
	"timer",
	/* local communication */
	"onNodeOnesided", "onNodePairwise",
	"onNodeOnesidedMinimum", "onNodePairwiseMinimum",
	/* remote communication */
	"offNodeOnesided", "offNodePairwise",
	"offNodeOnesidedMinimum", "offNodePairwiseMinimum"
};

/**
 \brief Create the measurement struct for the test
 \param tst Will tell the test how many times to run
 \param label A label for the measurement struct
*/
measurement_p net_measurement_create(test_p tst, char *label) {
	int i;
	measurement_p m = measurement_real_create(tst, label, NET_LEN);
	for (i = 0; i < NET_LEN; i++)
		strncpy(m->hist[i].label,net_labels[i],LABEL_LEN);
	return m;
}

/**
 \brief Exchanges messages between all the nodes on the system to test the network connections between them
 \param tst Tells how many cycles to run the test
 \param m Collects measurement data from the test
*/
void net_MPI_test(test_p tst, measurement_p m) {
#ifndef SHMEM
	gni_return_t status;
	gni_mem_handle_t sbuf_hndl, rbuf_hndl, rbuf_hndl_peer;
	gni_post_descriptor_t rma_desc;
	buffer_t *sbuf, *rbuf;
	double *cos, *cpw, *t;
	int i, icycle, istage, ierr, partner_rank;
	ORB_t t1, t2, t3;
	MPI_Status mpistatus;
	sbuf = comm_newbuffer(m->buflen);	/* exchange buffers */
	rbuf = comm_newbuffer(m->buflen);

	status = GNI_MemRegister(nic_handle,
				 (uint64_t)sbuf,
				 m->buflen,
				 NULL,
				 GNI_MEM_READWRITE,
				 -1,
				 &sbuf_hndl);
	assert(status == GNI_RC_SUCCESS);

	status = GNI_MemRegister(nic_handle,
				 (uint64_t)rbuf,
				 m->buflen,
				 NULL,
				 GNI_MEM_READWRITE,
				 -1,
				 &rbuf_hndl);
	assert(status == GNI_RC_SUCCESS);

	/*
	 * exchange memory handles for rbuf
	 */

	ierr = MPI_Sendrecv(&rbuf_hndl,
			    sizeof(rbuf_hndl),
			    MPI_BYTE,
			    partner_rank,
			    1010,
			    &rbuf_hndl_peer,
			    sizeof(rbuf_hndl),
			    MPI_BYTE,
			    partner_rank,
			    MPI_COMM_WORLD,
			    MPI_STATUS_IGNORE);

	/*
	 * code up the rma descriptor
	 */

	rma_desc.type = GNI_POST_RDMA_PUT;
	rma_desc.cq_mode = GNI_CQMODE_GLOBAL_EVENT |
				GNI_CQMODE_REMOTE_EVENT;
	rma_desc.dlvr_mode = dlvr_mode;
	rma_desc.local_addr = (uint64_t) send_buffer;
	rma_desc.local_mem_hndl = sbuf_hndl;
	rma_desc.remote_mem_hndl = rbuf_hndl_peer;
	rma_desc.length = m->buflen;
	rma_desc.src_cq_hndl = aft_nic.tx_cq;
	rma_desc.rdma_mode = 0;
	rma_desc.post_id = (uint64_t) &rma_desc;

	cos = (double *)malloc(tst->num_messages * sizeof(double));	/* array for onesided kernel timings */
	assert(cos != NULL);
	cpw = (double *)malloc(tst->num_messages * sizeof(double));	/* array for pairwise kernel timings */
	assert(cpw != NULL);
	t = (double *)malloc(tst->num_messages * sizeof(double));	/* array for timer overhead timings */
	assert(t != NULL);
	/* calibrate timer */
	ORB_calibrate();
	/* pre-synchronize all tasks */
	ierr = MPI_Barrier(MPI_COMM_WORLD);
	/*****************************************************************************
	 * A full set of samples for this task consists of message exchanges with each
	 * possible partner. The innermost loop below exchanges some number of messages
	 * between a particular pairing of partners. The middle loop steps through the
	 * possible partners. While the outmost allows us to aggregate multiple sets
	 * of samples to increase the total number of samples.
	 *****************************************************************************/
	for (icycle = 0; icycle < tst->num_cycles; icycle++) {
		/* step through the stage schedule */
		for (istage = 0; istage < tst->num_stages; istage++) {
			/* who's my buddy for this stage? */
			partner_rank = my_rank ^ istage;
			/* valid pairing */
			if ((partner_rank < num_ranks) && (partner_rank != my_rank)) {
				/* valid pair, proceed with test */
				ierr = 0;
				/***************************************/
				/* warm-up / pre-synchronize this pair */
				/***************************************/
				for (i = 0; i < tst->num_warmup; i++) {
					ORB_read(t1);
					ORB_read(t2);
					ierr += MPI_Sendrecv(sbuf->data, m->buflen, MPI_BYTE, partner_rank, 0,
							     rbuf->data, m->buflen, MPI_BYTE, partner_rank, 0,
							     MPI_COMM_WORLD, &mpistatus);
					ORB_read(t3);
				}
				assert(ierr == 0);
				/************************************************************/
				/* BEGIN PERFORMANCE KERNEL -- gather samples for this pair */
				/************************************************************/
				for (i = 0; i < tst->num_messages; i++) {
					/* for timer overhead estimate */
					ORB_read(t1);
					ORB_read(t2);
					/***************************************/
					/* begin timed communication primitive */
					/***************************************/

					status = GNI_PostRdma(aft_peer[partner_rank].ep,
		¬       ¬       ¬       ¬       ¬             &rma_desc);
					if (status != GNI_RC_SUCCESS) {
						fprintf(stderr,"SYSCONFIDENCE: (%s, PE %d) GNI_PostRdma returned %s\n",
							aft_my_cname, partner_rank, gni_err_str[status]);
					}

					/*
					 * wait for local side of transfer to complete
					 */
					status = GNI_RC_NOT_DONE;
					while (status != GNI_RC_SUCCESS)
						status = GNI_CqGetEvent(aft_nic.tx_cq,
									&current_event);
						check_tx_cqe(current_event);
					}

					/*
					 * wait for rx CQE from peer
					 */

					status = GNI_RC_NOT_DONE;
					while (status != GNI_RC_SUCCESS)
						status = GNI_CqGetEvent(aft_nic.rx_cq,
									&current_event);
						check_rx_cqe(current_event);
					}

					/*************************************/
					/* end timed communication primitive */
					/*************************************/
					ORB_read(t3);
					/* save the timings */
					t[i] = ORB_seconds(t2, t1);
					cos[i] = ORB_seconds(t3, t2);
				}
				/************************************************************/
				/* END PERFORMANCE KERNEL -- samples gathered for this pair */
				/************************************************************/
				assert(ierr == 0);
				/* exchange array of local timings with partner */
				ierr += MPI_Sendrecv(cos, tst->num_messages, MPI_DOUBLE, partner_rank, 0,
						     cpw, tst->num_messages, MPI_DOUBLE, partner_rank, 0,
						     MPI_COMM_WORLD, &mpistatus);
				assert(ierr == 0);
				/* pairwise as average, comparable to one-sided */
				for (i = 0; i < tst->num_messages; i++) {
					cpw[i] = (cpw[i] + cos[i]) / 2.0;
				}
				/* bin the t, cos, and cpw results for this cycle of this pair in p */
				net_measurement_bin(tst, m, t, cos, cpw, (node_id[my_rank] == node_id[partner_rank]));
			} /* if valid pairing */
		} /* for istage */
	} /* for icycle */
	free(t);
	free(cpw);
	free(cos);
	comm_freebuffer(rbuf);
	comm_freebuffer(sbuf);
	status = GNI_MemDeregister(aft_nic.nic, sbuf_hndl);
	assert(status == GNI_RC_SUCCESS);
	status = GNI_MemDeregister(aft_nic.nic, rbuf_hndl);
	assert(status == GNI_RC_SUCCESS);

#endif
	return;
}

/**
 \brief Converts the time measurements to bin
*/
void net_measurement_bin(test_p tst, measurement_p m, double *t, double *cos, double *cpw, int LOCAL) {
	double cosmin, cpwmin;
	uint64_t *ptimd, *posd, *ppwd, *posmd, *ppwmd;
	int i;
	if (LOCAL) { /* bin these values as local communication */
		posd = m->hist[onNodeOnesided].dist;
		ppwd = m->hist[onNodePairwise].dist;
		posmd = m->hist[onNodeOnesidedMinimum].dist;
		ppwmd = m->hist[onNodePairwiseMinimum].dist;
	} else { /* bin these values as remote communication */
		posd = m->hist[offNodeOnesided].dist;
		ppwd = m->hist[offNodePairwise].dist;
		posmd = m->hist[offNodeOnesidedMinimum].dist;
		ppwmd = m->hist[offNodePairwiseMinimum].dist;
	}
	ptimd = m->hist[timer].dist;
	cosmin = cpwmin = 1.0e+16;
	for (i = 0; i < tst->num_messages; i++) {
		/* bin the individual results */
		if (t != NULL) {
			if (t[i] >= 0.0)
				ptimd[time2bin(tst,t[i])]++;
		}
		if (cos[i] >= 0.0)
			posd[time2bin(tst,cos[i])]++;
		if (cpw[i] >= 0.0)
			ppwd[time2bin(tst,cpw[i])]++;
		/* save the minimums for now */
		if ((cos[i] > 0.0) && (cos[i] < cosmin))
			cosmin = cos[i];
		if ((cpw[i] > 0.0) && (cpw[i] < cpwmin))
			cpwmin = cpw[i];
	}
	/* now bin the minimums for this communications pair */
	if (cosmin > 0.0)
		posmd[time2bin(tst,cosmin)]++;
	if (cpwmin > 0.0)
		ppwmd[time2bin(tst,cpwmin)]++;
}


