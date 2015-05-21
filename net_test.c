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

#ifdef SHMEM
	#include <mpp/shmem.h>
#else
	#include <mpi.h>
#endif

#ifdef UGNI_DIRECT
#include "ugni_utils.h"
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
void net_SHMEM_test(test_p tst, measurement_p m) {
#ifdef SHMEM
	static int sync;
	sync = my_rank;
	buffer_t *sbuf, *rbuf;
	double *cos, *cpw, *t;
	int i, icycle, istage, partner_rank;
	ORB_t t1, t2, t3;
//Make
	sbuf = comm_newbuffer(m->buflen);	/* exchange buffers */
	rbuf = comm_newbuffer(m->buflen);
	cos = (double *)shmalloc(tst->num_messages * sizeof(double));	/* array for onesided kernel timings */
	assert(cos != NULL);
	cpw = (double *)malloc(tst->num_messages * sizeof(double));	/* array for pairwise kernel timings */
	assert(cpw != NULL);
	t = (double *)malloc(tst->num_messages * sizeof(double));	/* array for timer overhead timings */
	assert(t != NULL);

// Exec
	/* calibrate timer */
	ORB_calibrate();
	/* pre-synchronize all tasks */
	shmem_barrier_all();
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
			shmem_barrier_all();
			/* who's my buddy for this stage? */
			partner_rank = my_rank ^ istage;
			/* valid pairing */
			if ((partner_rank < num_ranks) && (partner_rank != my_rank)) {
				/* valid pair, proceed with test */
				
				/***************************************/
				/* warm-up / pre-synchronize this pair */
				/***************************************/
				for (i = 0; i < tst->num_warmup; i++) {
					ORB_read(t1);
					ORB_read(t2);
					shmem_getmem(rbuf->data, sbuf->data, m->buflen, partner_rank);
					ORB_read(t3);
				}
				
				/* synchronize partners */
				shmem_int_p(&sync, my_rank, partner_rank);
				shmem_int_wait_until(&sync, SHMEM_CMP_EQ, partner_rank);
				sync = my_rank;
				
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
					shmem_getmem(rbuf->data, sbuf->data, m->buflen, partner_rank);
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
				
				/* ensure partner has completed sample collection */
				shmem_int_p(&sync, my_rank, partner_rank);
				shmem_int_wait_until(&sync, SHMEM_CMP_EQ, partner_rank);
				sync = my_rank;
				
				/* get partner's array of local timings */
				shmem_double_get(cpw, cos, tst->num_messages, partner_rank);
				
				/* pairwise as average, comparable to one-sided */
				for (i = 0; i < tst->num_messages; i++) {
					cpw[i] = (cpw[i] + cos[i]) / 2.0;
				}
				/* bin the t, cos, and cpw results for this cycle of this pair in p */
				net_measurement_bin(tst, m, t, cos, cpw, (node_id[my_rank] == node_id[partner_rank]));
			} /* if valid pairing */
		} /* for istage */
	} /* for icycle */

// Kill
	free(t);
	free(cpw);
	shfree(cos);
	comm_freebuffer(rbuf);
	comm_freebuffer(sbuf);
#endif
	return;
}

/**
 \brief Exchanges messages between all the nodes on the system to test the network connections between them
 \param tst Tells how many cycles to run the test
 \param m Collects measurement data from the test
*/
void net_MPI_test(test_p tst, measurement_p m) {
#ifndef SHMEM
	buffer_t *sbuf, *rbuf;
	double *cos, *cpw, *t;
	int i, icycle, istage, ierr, partner_rank;
	ORB_t t1, t2, t3;
	MPI_Status mpistatus;
	sbuf = comm_newbuffer(m->buflen);	/* exchange buffers */
	rbuf = comm_newbuffer(m->buflen);
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
					ierr += MPI_Sendrecv(sbuf->data, m->buflen, MPI_BYTE, partner_rank, 0,
							     rbuf->data, m->buflen, MPI_BYTE, partner_rank, 0,
							     MPI_COMM_WORLD, &mpistatus);
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
#endif
	return;
}

/**
 \brief Exchanges messages between all the nodes on the system to test the network connections between them
 \param tst Tells how many cycles to run the test
 \param m Collects measurement data from the test
*/
void net_uGNI_test(test_p tst, measurement_p m) {
#ifdef UGNI_DIRECT
	gni_return_t status;
	gni_mem_handle_t sbuf_hndl, rbuf_hndl, rbuf_hndl_peer;
	gni_cq_entry_t current_event;
	gni_post_descriptor_t rma_desc,*rma_desc_ptr;
	buffer_t *sbuf, *rbuf;
	double *cos, *cpw, *t;
	int i, icycle, istage, ierr, partner_rank;
	ORB_t t1, t2, t3;
	MPI_Status mpistatus;
	gni_ep_handle_t gni_ep;
	typedef struct {
		gni_mem_handle_t mdh;
		uint64_t	addr;
	} mdh_addr_t;

	mdh_addr_t my_info, partner_info;

	sbuf = comm_newbuffer(m->buflen);	/* exchange buffers */
	rbuf = comm_newbuffer(m->buflen);

#if 1
	fprintf(stderr,"ENTERING uGNI TEST %d\n", m->buflen);
#endif
	/*
	 * create and endpoint and
	 * bind to peer
	 */

	status = GNI_EpCreate(gni_nic.nic,
				gni_nic.tx_cq,
				&gni_ep);
	assert(status == GNI_RC_SUCCESS);


	status = GNI_MemRegister(gni_nic.nic,
				 (uint64_t)sbuf->data,
				 m->buflen,
				 NULL,
				 GNI_MEM_READWRITE,
				 -1,
				 &sbuf_hndl);
	assert(status == GNI_RC_SUCCESS);

	status = GNI_MemRegister(gni_nic.nic,
				 (uint64_t)rbuf->data,
				 m->buflen,
				 gni_nic.rx_cq,
				 GNI_MEM_READWRITE,
				 -1,
				 &rbuf_hndl);
	assert(status == GNI_RC_SUCCESS);

	/*
	 * code up the rma descriptor
	 */

	rma_desc.type = GNI_POST_RDMA_PUT;
	rma_desc.cq_mode = GNI_CQMODE_GLOBAL_EVENT |
				GNI_CQMODE_REMOTE_EVENT;
	rma_desc.dlvr_mode = sysconf_ugni_dlvr_mode;
	rma_desc.local_addr = (uint64_t) sbuf->data;
	rma_desc.local_mem_hndl = sbuf_hndl;
	rma_desc.length = m->buflen;
	rma_desc.src_cq_hndl = gni_nic.tx_cq;
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

			my_info.mdh = rbuf_hndl;
			my_info.addr = (uint64_t)rbuf->data;
#if 0
			fprintf(stderr,"rank %d settting info vaddr 0x%016lx 0x%016lx 0x%016lx\n",
				my_rank, (uint64_t)rbuf->data, my_info.mdh.qword1, my_info.mdh.qword2);
#endif

			ierr = MPI_Sendrecv(&my_info,
				    sizeof(mdh_addr_t),
				    MPI_BYTE,
				    partner_rank,
				    1010,
				    &partner_info,
				    sizeof(mdh_addr_t),
				    MPI_BYTE, partner_rank,
				    1010,
				    MPI_COMM_WORLD,
				    &mpistatus);

#if 0
			fprintf(stderr,"rank %d has partner vaddr 0x%016lx 0x%016lx 0x%016lx\n",
				my_rank, (uint64_t)partner_info.addr, partner_info.mdh.qword1,
				partner_info.mdh.qword2);
#endif

			rma_desc.remote_addr = (uint64_t) partner_info.addr;
			rma_desc.remote_mem_hndl = partner_info.mdh;

			status = GNI_EpBind(gni_ep,
					sysconf_peer_gni_info[partner_rank].gni_addr,
					partner_rank);
#if 0
			fprintf(stderr,"GNI_EPBind returned %s\n",gni_err_str[status]);
#endif
			assert(status == GNI_RC_SUCCESS);

			/* valid pairing */
			if ((partner_rank < num_ranks) && (partner_rank !=
								my_rank)) {
				/* valid pair, proceed with test */
				ierr = 0;
				/***************************************/
				/* warm-up / pre-synchronize this pair */
				/***************************************/
				for (i = 0; i < tst->num_warmup; i++) {
					ORB_read(t1);
					ORB_read(t2);

					status = GNI_EpSetEventData(gni_ep,
								    my_rank,
								    partner_rank + i);

					status = GNI_PostRdma(gni_ep,
								&rma_desc);
					if (status != GNI_RC_SUCCESS) {
						fprintf(stderr,"SYSCONFIDENCE: (%s, PE %d) GNI_PostRdma returned %s\n",
							sysconf_my_gni_info.cname,
							partner_rank, gni_err_str[status]);
					}

#if 0
					fprintf(stderr,"rank %d Waiting for local completion\n", my_rank);
#endif
					/*
					 * wait for local side of transfer to complete
					 */
					status = GNI_RC_NOT_DONE;
					while (status != GNI_RC_SUCCESS) {
						status = GNI_CqGetEvent(gni_nic.tx_cq,
									&current_event);
						if (status == GNI_RC_TRANSACTION_ERROR) {
						    fprintf(stderr,"Got a network cqe error waiting for put from %s to %s\n",
							sysconf_my_gni_info.cname,
							sysconf_peer_gni_info[partner_rank].cname);
						}
					}

					assert(GNI_CQ_GET_TYPE(current_event) ==
						GNI_CQ_EVENT_TYPE_POST);

					status = GNI_GetCompleted(gni_nic.tx_cq,
								  current_event,
								  &rma_desc_ptr);
					assert((status == GNI_RC_SUCCESS) ||
					       (status == GNI_RC_TRANSACTION_ERROR));

					if (status == GNI_RC_TRANSACTION_ERROR){
						fprintf(stderr,"Got a network error putting to %s from %s\n",
							sysconf_peer_gni_info[partner_rank].cname,
							sysconf_my_gni_info.cname);
					}

#if 0
					fprintf(stderr,"rank %d Waiting for remote completion\n",my_rank);
#endif
					/*
					 * wait for rx CQE from peer
					 */

					status = GNI_RC_NOT_DONE;
					while (status != GNI_RC_SUCCESS) {
						status = GNI_CqGetEvent(gni_nic.rx_cq,
									&current_event);
						if (status == GNI_RC_TRANSACTION_ERROR) {
						    fprintf(stderr,"Got a network cqe error waiting for put from %s to %s\n",
							sysconf_my_gni_info.cname,
							sysconf_peer_gni_info[partner_rank].cname);
						}
					}

					ORB_read(t3);
				}

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

					status = GNI_EpSetEventData(gni_ep,
								    my_rank,
								    partner_rank + i);

					status = GNI_PostRdma(gni_ep,
								&rma_desc);
					if (status != GNI_RC_SUCCESS) {
						fprintf(stderr,"SYSCONFIDENCE: (%s, PE %d) GNI_PostRdma returned %s\n",
							sysconf_my_gni_info.cname, partner_rank, gni_err_str[status]);
					}

					/*
					 * wait for local side of transfer to complete
					 */
#if 0
					fprintf(stderr,"rank %d Waiting for local completion\n", my_rank);
#endif
					status = GNI_RC_NOT_DONE;
					while (status != GNI_RC_SUCCESS) {
						status = GNI_CqGetEvent(gni_nic.tx_cq,
									&current_event);
						if (status == GNI_RC_TRANSACTION_ERROR) {
						    fprintf(stderr,"Got a network cqe error waiting for put from %s to %s\n",
							sysconf_my_gni_info.cname,
							sysconf_peer_gni_info[partner_rank].cname);
						}
					}

					assert(GNI_CQ_GET_TYPE(current_event) ==
						GNI_CQ_EVENT_TYPE_POST);

					status = GNI_GetCompleted(gni_nic.tx_cq,
								  current_event,
								  &rma_desc_ptr);
					assert((status == GNI_RC_SUCCESS) ||
					       (status == GNI_RC_TRANSACTION_ERROR));

					if (status == GNI_RC_TRANSACTION_ERROR){
						fprintf(stderr,"Got a network error putting to %s from %s\n",
							sysconf_peer_gni_info[partner_rank].cname,
							sysconf_my_gni_info.cname);
					}


					/*
					 * wait for rx CQE from peer
					 */

#if 0
					fprintf(stderr,"rank %d Waiting for remote completion\n",my_rank);
#endif
					status = GNI_RC_NOT_DONE;
					while (status != GNI_RC_SUCCESS) {
						status = GNI_CqGetEvent(gni_nic.rx_cq,
									&current_event);
						if (status == GNI_RC_TRANSACTION_ERROR) {
						    fprintf(stderr,"Got a network cqe error waiting for put from %s to %s\n",
							sysconf_my_gni_info.cname,
							sysconf_peer_gni_info[partner_rank].cname);
						}
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
			status = GNI_EpUnbind(gni_ep);
			assert(status == GNI_RC_SUCCESS);
		} /* for istage */
	} /* for icycle */
	free(t);
	free(cpw);
	free(cos);
	comm_freebuffer(rbuf);
	comm_freebuffer(sbuf);
	status = GNI_MemDeregister(gni_nic.nic, &sbuf_hndl);
	assert(status == GNI_RC_SUCCESS);
	status = GNI_MemDeregister(gni_nic.nic, &rbuf_hndl);
	assert(status == GNI_RC_SUCCESS);

	status = GNI_EpDestroy(gni_ep);
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


