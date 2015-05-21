
#include "aft_internal.h"

aft_nic_t aft_nic;
aft_conn_info_t my_conn_info;
gni_ep_handle_t *aft_ep_hndls;

/*
 * by default we use in order routing on aries since
 * first order of business for sysconfidence test
 * is to very link status - which is easier when
 * using in order.
 */

uint16_t sysconf_gni_dlvr_mode = GNI_DLVMODE_IN_ORDER;

/*
 * we ask for fma sharing - hopefully won't stress
 * out craypich
 */
uint32_t sysconf_gni_cdm_modes = GNI_CDM_MODE_FMA_SHARED;


int ugni_init(int cdm_modes)
{
	char *env_str;
	int i, rc, my_rank, nranks;
	int device_id = 0; /* only 1 aries nic/node */
	uint32_t bytes_per_smsg;
	gni_smsg_attr_t smsg_attr;
	gni_smsg_attr_t *all_smsg_attrs = NULL;

	/*
	 * we need pmi to be fired up in order to initialize
	 */
	rc = PMI_Get_rank(&my_rank);
	if (rc != PMI_SUCCESS) {
		fprintf(stderr, "SYSCONF: PMI_Get_rank failed %d\n", my_rank);
		exit(-1); /* exit here since likely MPI_Abort isn't going to work either */
	}

	/*
	 * parse the environment
	 */

	/*
	 * see gni_pub.h
	 */
	env_str = getenv("SYSCONF_DLVR_MODE");
	if (env_str != NULL) {
		if (strcmp(env_str,"GNI_DLVMODE_NMIN_HASH")) {
			sysconf_gni_dlvr = GNI_DLVMODE_NMIN_HASH;
		}
		else if (strcmp(value,"GNI_DLVMODE_MIN_HASH")) {
			sysconf_gni_dlvr = GNI_DLVMODE_MIN_HASH;
		}
		else if (strcmp(value,"GNI_DLVMODE_ADAPTIVE0")) {
			sysconf_gni_dlvr = GNI_DLVMODE_ADAPTIVE0;
		}
		else if (strcmp(value,"GNI_DLVMODE_ADAPTIVE1")) {
			sysconf_gni_dlvr = GNI_DLVMODE_ADAPTIVE1;
		}
		else if (strcmp(value,"GNI_DLVMODE_ADAPTIVE2")) {
			sysconf_gni_dlvr = GNI_DLVMODE_ADAPTIVE2;
		}
		else if (strcmp(value,"GNI_DLVMODE_ADAPTIVE3")) {
			sysconf_gni_dlvr = GNI_DLVMODE_ADAPTIVE3;
		}
        }

	if (getenv("SYSCONF_USE_PCI_IOMMU") != NULL) {
		sysconf_ugni_cdm_modes |= GNI_CDM_MODE_USE_PCI_IOMMU;
	}

	if (getenv("SYSCONF_FLBTE_DISABLE") != NULL) {
		sysconf_ugni_cdm_modes |= GNI_CDM_MODE_FLBTE_DISABLE;
	}

	if (getenv("SYSCONF_BTE_SINGLE_CHANNEL") != NULL) {
		sysconf_ugni_cdm_modes |= GNI_CDM_MODE_BTE_SINGLE_CHANNEL;
	}

	if (getenv("SYSCONF_UGNI_DISPLAY") != NULL) && (my_rank == 0)) {
		fprintf(stderr,"SYSCONF: UGNI DLVR MODE - ");
		switch (sysconf_gni_dlvr_mode) {
		GNI_DLVMODE_IN_ORDER:
			fprintf(stderr,"GNI_DLVMODE_IN_ORDER\n");
			break;
		GNI_DLVMODE_NMIN_HASH:
			fprintf(stderr,"GNI_DLVMODE_NMIN_HASH\n");
			break;
		GNI_DLVMODE_MIN_HASH:
			fprintf(stderr,"GNI_DLVMODE_MIN_HASH\n");
			break;
		GNI_DLVMODE_ADAPTIVE0:
			fprintf(stderr,"GNI_DLVMODE_ADAPTIVE0\n");
			break;
		GNI_DLVMODE_ADAPTIVE1:
			fprintf(stderr,"GNI_DLVMODE_ADAPTIVE1\n");
			break;
		GNI_DLVMODE_ADAPTIVE2:
			fprintf(stderr,"GNI_DLVMODE_ADAPTIVE2\n");
			break;
		GNI_DLVMODE_ADAPTIVE3:
			fprintf(stderr,"GNI_DLVMODE_ADAPTIVE3\n");
			break;
		}
		if (sysconf_ugni_cdm_modes & GNI_CDM_MODE_USE_PCI_IOMMU) {
			fprintf(stderr,"SYSCONF: USING PCI IOMMU\n");
		}
		if (sysconf_ugni_cdm_modes & GNI_CDM_MODE_FLBTE_DISABLE) {
			fprintf(stderr,"SYSCONF: FLBTE disabled\n");
		}
		if (sysconf_ugni_cdm_modes & GNI_CDM_MODE_BTE_SINGLE_CHANNEL) {
			fprintf(stderr,"SYSCONF: Using BTE in single channel mode\n");
		}
	}

	/*
	 * gather up cnames
	 */

	rc = PMI_Allgather(&my_conn_info,
			   all_smsg_attrs,
                           sizeof(my_conn_info));
	if (rc != PMI_SUCCESS) {
		fprintf(stderr,"SYSCONF: PMI_Allgather returned %d aborting...\n",rc);
		MPI_Abort(MPI_COMM_WORLD,-1);
	}
	rc = PMI_Get_size(&nranks);
	if (rc != PMI_SUCCESS) {
		AFT_WARN("PMI_Get_size returned %d\n",rc);
		goto err;
	}

	/*
	 * Get the GNI RDMA credentials from PMI - we use the first
	 * ptag/cookie since at some point one might want to use the
	 * shmem and ugni together
	 */

	ptag = __get_ptag();
	cookie = __get_cookie();

	status = GNI_CdmCreate(my_rank,
			       ptag,
			       cookie,
			       cdm_modes,
			       &aft_nic.cdm);
	if (status != GNI_RC_SUCCESS) {
		AFT_WARN("GNI_CdmCreate returned %s\n",
			gni_err_str[status]);
		goto err;
	}

	status = GNI_CdmAttach(aft_nic.cdm,
			       device_id,
			       &local_address,
			       &aft_nic.nic);
	if (status != GNI_RC_SUCCESS) {
		AFT_WARN("GNI_CdmAttach returned %s\n",
			gni_err_str[status]);
		goto err1;
	}

	/*
	 * create a TX CQ
	 */

	status = GNI_CqCreate(aft_nic.nic,
			      number_of_cq_entries,
			      0,
			      GNI_CQ_NOBLOCK | GNI_CQ_PHYS_PAGES,
			      NULL,
			      NULL,
			      &aft_nic.tx_cq);
	if (status != GNI_RC_SUCCESS) {
		AFT_WARN("GNI_CqCreate returned %s\n",
			gni_err_str[status]);
		goto err1;
	}

	/*
	 * create a RX CQ
	 */
        status = GNI_CqCreate(aft_nic.nic,
			      number_of_rx_cq_entries,
			      GNI_CQ_NOBLOCK | GNI_CQ_PHYS_PAGES,
			      NULL,
			      NULL,
			      &aft_nic.rx_cq);
	if (status != GNI_RC_SUCCESS) {
		AFT_WARN("GNI_CqCreate returned %s\n",
			gni_err_str[status]);
		goto err2;
	}

	rc = PMI_Allgather(&my_conn_info,
			   all_smsg_attrs,
                           sizeof(my_conn_info));
	if (rc != PMI_SUCCESS) {
		AFT_WARN("PMI_Allgather returned %d\n", rc);
		goto err3;
	}

	/*
	 * Set up the endpoints
	 */

	aft_ep_hndls  = malloc(nranks * sizeof(gni_ep_handle_t));
	if (ep_hndls == NULL) {
		AFT_WARN("calloc of ep_hndls failed\n");
		goto err3;
	}

	i = 0;
	while(i < nranks) {
		the_rank = all_smsg_attrs[i].my_rank;
        	status = GNI_EpCreate(aft_nic.nic,
					aft_nic.tx_cq,
					&aft_ep_hndls[the_rank]);
		if (status != GNI_RC_SUCCESS) {
			AFT_WARN("GNI_MemRegister returned %s\n",
				gni_err_str[status]);
			goto err3;
		}
		all_smsg_attrs[i].attr.mbox_offset = bytes_per_smsg * the_rank;

		smsg_attr.msg_type = GNI_SMSG_TYPE_MBOX_AUTO_RETRANSMIT;
		smsg_attr.mbox_maxcredit = 16;
		smsg_attr.msg_maxsize = 512;
		smsg_attr.mbox_offset = bytes_per_smsg * the_rank;

		status = GNI_SmsgInit(aft_ep_hndls[the_rank],
					&smsg_attr,
					&all_smsg_attrs[i]);
		if (status != GNI_RC_SUCCESS) {
			AFT_WARN("GNI_MemRegister returned %s\n",
				gni_err_str[status]);
			goto err3;
		}
		++i;;
	}

	/*
	 * need to barrier here to make sure all ranks have
	 * initialized their endpoints for messaging
	 */

	rc = PMI_Barrier();
	if (rc != PMI_SUCCESS) {
		AFT_WARN("GNI_MemRegister returned %s\n",
			gni_err_str[status]);
		goto err;
	}

	return 0;

err3:
	GNI_CqDestroy(aft_nic.rx_cq);
err2:
	GNI_CqDestroy(aft_nic.tx_cq);
err1:
	GNI_CdmDestroy(aft_nic.cdm);
err:
	if (all_smsg_attrs != NULL)
		free(all_smsg_attrs);
	if (aft_ep_hndls != NULL)
		free(aft_ep_hndls);

	PMI_Finalize();
	return -1;
}

int aft_finalize(int cdm_modes)
{
	int i, rc, my_rank, nranks;

	/*
	 * wasn't initialized, nothing to do
	 */

	if (!aft_initialized)
		return -1;

	/*
	 * PMI barrier to make sure everyone's done with
	 * remote memory access
	 */

	rc = PMI_Barrier();
	if (rc != PMI_SUCCESS) {
		AFT_WARN("PMI_Init returned %d\n",rc);
		goto err;
	}

	/*
	 * destroy the endpoints
	 */

	for (i = 0; i < nranks; i++) {
		status = GNI_EpDestroy(aft_ep_hndls[i]);
		if (status != GNI_RC_SUCCESS) {
			AFT_WARN("GNI_EpDestroy returned %s\n",
				gni_err_str[status]);
			goto err;
		}
	}

	status = GNI_CqDestroy(aft_nic.tx_cq);
	if (status != GNI_RC_SUCCESS) {
		AFT_WARN("GNI_CqDestroy returned %s\n",
			gni_err_str[status]);
		goto err;
	}

	status = GNI_CqDestroy(aft_nic.rx_cq);
	if (status != GNI_RC_SUCCESS) {
		AFT_WARN("GNI_CqDestroy returned %s\n",
			gni_err_str[status]);
		goto err;
	}

	status = GNI_CdmDestroy(aft_nic.cdm);
	if (status != GNI_RC_SUCCESS) {
		AFT_WARN("GNI_CqDestroy returned %s\n",
			gni_err_str[status]);
		goto err;
	}

	rc = PMI_Finalize();
	if (rc != PMI_SUCCESS)
		AFT_WARN("PMI_Finalize returned %d\n", rc);

	return 0;

err:
	return -1;
}
