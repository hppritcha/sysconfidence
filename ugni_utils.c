
#include "ugni_utils.h"
#include "mpi.h"
#include <errno.h>

sysconf_gni_nic_t gni_nic;
sysconf_gni_peer_info_t *sysconf_peer_gni_info;
sysconf_gni_peer_info_t sysconf_my_gni_info;

static int ugni_initialized;
/*
 * by default we use in order routing on aries since
 * first order of business for sysconfidence test
 * is to very link status - which is easier when
 * using in order.
 */

uint16_t sysconf_ugni_dlvr_mode = GNI_DLVMODE_IN_ORDER;

/*
 * we ask for fma sharing - hopefully won't stress
 * out craypich
 */
static uint32_t ugni_cdm_modes = GNI_CDM_MODE_FMA_SHARED;

/*
 * forward prototypes
 */

/**
 * @brief get local cname of node from /proc/cray_xt
 * @param[out] cname  pointer to buffer where cname will be returned
 * @param[in] len     length of buffer where cname will be returned
 *
 * @return 0 on success, -1 on failure
 */
static int _get_cname(char *cname, int len);

static int _get_ptag(uint8_t *ptag);
static int _get_cookie(uint32_t *cookie);

int ugni_init(void)
{
	gni_return_t status;
	char *env_str;
	uint8_t ptag;
	uint32_t cookie;
	int i, rc, my_rank, nranks;
	int device_id = 0; /* only 1 aries nic/node */
	sysconf_gni_peer_info_t my_gni_info;
	sysconf_gni_peer_info_t *tmp_info;

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
			sysconf_ugni_dlvr_mode = GNI_DLVMODE_NMIN_HASH;
		} else if (strcmp(env_str,"GNI_DLVMODE_MIN_HASH")) {
			sysconf_ugni_dlvr_mode = GNI_DLVMODE_MIN_HASH;
		} else if (strcmp(env_str,"GNI_DLVMODE_ADAPTIVE0")) {
			sysconf_ugni_dlvr_mode = GNI_DLVMODE_ADAPTIVE0;
		} else if (strcmp(env_str,"GNI_DLVMODE_ADAPTIVE1")) {
			sysconf_ugni_dlvr_mode = GNI_DLVMODE_ADAPTIVE1;
		}else if (strcmp(env_str,"GNI_DLVMODE_ADAPTIVE2")) {
			sysconf_ugni_dlvr_mode = GNI_DLVMODE_ADAPTIVE2;
		} else if (strcmp(env_str,"GNI_DLVMODE_ADAPTIVE3")) {
			sysconf_ugni_dlvr_mode = GNI_DLVMODE_ADAPTIVE3;
		}
        }

	if (getenv("SYSCONF_USE_PCI_IOMMU") != NULL) {
		ugni_cdm_modes |= GNI_CDM_MODE_USE_PCI_IOMMU;
	}

	if (getenv("SYSCONF_FLBTE_DISABLE") != NULL) {
		ugni_cdm_modes |= GNI_CDM_MODE_FLBTE_DISABLE;
	}

	if (getenv("SYSCONF_BTE_SINGLE_CHANNEL") != NULL) {
		ugni_cdm_modes |= GNI_CDM_MODE_BTE_SINGLE_CHANNEL;
	}

	if ((getenv("SYSCONF_UGNI_DISPLAY") != NULL) && (my_rank == 0)) {
		fprintf(stderr,"SYSCONF: UGNI DLVR MODE - ");
		if (sysconf_ugni_dlvr_mode == GNI_DLVMODE_IN_ORDER) {
			fprintf(stderr,"GNI_DLVMODE_IN_ORDER\n");
		} else if (sysconf_ugni_dlvr_mode == GNI_DLVMODE_NMIN_HASH) {
			fprintf(stderr,"GNI_DLVMODE_NMIN_HASH\n");
		} else if (sysconf_ugni_dlvr_mode == GNI_DLVMODE_MIN_HASH) {
			fprintf(stderr,"GNI_DLVMODE_MIN_HASH\n");
		} else if (sysconf_ugni_dlvr_mode == GNI_DLVMODE_ADAPTIVE0) {
			fprintf(stderr,"GNI_DLVMODE_ADAPTIVE0\n");
		} else if (sysconf_ugni_dlvr_mode == GNI_DLVMODE_ADAPTIVE1) {
			fprintf(stderr,"GNI_DLVMODE_ADAPTIVE1\n");
		} else if (sysconf_ugni_dlvr_mode == GNI_DLVMODE_ADAPTIVE2) {
			fprintf(stderr,"GNI_DLVMODE_ADAPTIVE2\n");
		} else if (sysconf_ugni_dlvr_mode == GNI_DLVMODE_ADAPTIVE3) {
			fprintf(stderr,"GNI_DLVMODE_ADAPTIVE3\n");
		}
		if (ugni_cdm_modes & GNI_CDM_MODE_USE_PCI_IOMMU) {
			fprintf(stderr,"SYSCONF: USING PCI IOMMU\n");
		}
		if (ugni_cdm_modes & GNI_CDM_MODE_FLBTE_DISABLE) {
			fprintf(stderr,"SYSCONF: FLBTE disabled\n");
		}
		if (ugni_cdm_modes & GNI_CDM_MODE_BTE_SINGLE_CHANNEL) {
			fprintf(stderr,"SYSCONF: Using BTE in single channel mode\n");
		}
	}

	sysconf_my_gni_info.rank = my_rank;

	rc = _get_cname(sysconf_my_gni_info.cname,
			sizeof(sysconf_my_gni_info.cname));
	if (rc != 0) {
		fprintf(stderr,"SYSCONF(%d): _get_cname failed\n",my_rank);
		MPI_Abort(MPI_COMM_WORLD,-1);
	}

	/*
 	 * get rdma credentials needed to initialize ugni cdm
 	 * use the same credentials as craypich/openmpi in case
 	 * wanting to also use shmem in the sysconfidence test
 	 * at the same time
 	 */

	rc = _get_ptag(&ptag);
	if (rc != 0) {
		fprintf(stderr,"SYSCONF(%s,%d): _get_ptag failed\n",
			sysconf_my_gni_info.cname,
			sysconf_my_gni_info.rank);
		MPI_Abort(MPI_COMM_WORLD,-1);
	}

	rc = _get_cookie(&cookie);
	if (rc != 0) {
		fprintf(stderr,"SYSCONF(%s,%d): _get_cookie failed\n",
			sysconf_my_gni_info.cname,
			sysconf_my_gni_info.rank);
		MPI_Abort(MPI_COMM_WORLD,-1);
	}

	/*
 	 * create the cdm
 	 */
	status = GNI_CdmCreate(sysconf_my_gni_info.rank,
				ptag,
				cookie,
				ugni_cdm_modes,
				&gni_nic.cdm);
	if (status != GNI_RC_SUCCESS) {
		fprintf(stderr,"SYSCONF(%s,%d): GNI_CdmCreate failed %s\n",
			sysconf_my_gni_info.cname,
			sysconf_my_gni_info.rank,
			gni_err_str[status]);
		MPI_Abort(MPI_COMM_WORLD,-1);
	}

	/*
 	 * attach cdm to the local nic
 	 */
	status = GNI_CdmAttach(gni_nic.cdm,
			       device_id,
			       &sysconf_my_gni_info.gni_addr,
			       &gni_nic.nic);
	if (status != GNI_RC_SUCCESS) {
		fprintf(stderr,"SYSCONF(%s,%d): GNI_CdmAttach failed %s\n",
			sysconf_my_gni_info.cname,
			sysconf_my_gni_info.rank,
			gni_err_str[status]);
		MPI_Abort(MPI_COMM_WORLD,-1);
	}

	/*
	 * create a TX CQ
	 */

	status = GNI_CqCreate(gni_nic.nic,
			      1024, /* plenty, we're doing simple comm. */
			      0,
			      GNI_CQ_NOBLOCK | GNI_CQ_PHYS_PAGES,
			      NULL,
			      NULL,
			      &gni_nic.tx_cq);
	if (status != GNI_RC_SUCCESS) {
		fprintf(stderr,"SYSCONF(%s,%d): GNI_CqCreate failed %s\n",
			sysconf_my_gni_info.cname,
			sysconf_my_gni_info.rank,
			gni_err_str[status]);
		MPI_Abort(MPI_COMM_WORLD,-1);
	}

	/*
	 * create a RX CQ
	 */
        status = GNI_CqCreate(gni_nic.nic,
			      10 * 1024, /* just in case we retransmit */
			      0,
			      GNI_CQ_NOBLOCK | GNI_CQ_PHYS_PAGES,
			      NULL,
			      NULL,
			      &gni_nic.rx_cq);
	if (status != GNI_RC_SUCCESS) {
		fprintf(stderr,"SYSCONF(%s,%d): GNI_CqCreate failed %s\n",
			sysconf_my_gni_info.cname,
			sysconf_my_gni_info.rank,
			gni_err_str[status]);
		MPI_Abort(MPI_COMM_WORLD,-1);
	}

	/*
 	 * exchange info
 	 */

	rc = PMI_Get_size(&nranks);
	if (rc != PMI_SUCCESS) {
		fprintf(stderr,"SYSCONF(%s,%d): PMI_Get_size failed %d\n",
			sysconf_my_gni_info.cname,
			sysconf_my_gni_info.rank);
		MPI_Abort(MPI_COMM_WORLD,-1);
	}

	tmp_info = malloc(sizeof(sysconf_gni_peer_info_t) * nranks);
	if (tmp_info == NULL) {
		fprintf(stderr,"SYSCONF(%s,%d): malloc of %d bytes failed\n",
			sysconf_my_gni_info.cname,
			sysconf_my_gni_info.rank,
			sizeof(sysconf_gni_peer_info_t) * nranks);
		MPI_Abort(MPI_COMM_WORLD,-1);
	}

	rc = PMI_Allgather(&sysconf_my_gni_info,
			   tmp_info,
                           sizeof(sysconf_gni_peer_info_t));
	/*
 	 * now rearrange the info into the final array - PMI_Allgather
 	 * doesn't necessarily deliver the data in rank order
 	 */

	sysconf_peer_gni_info = malloc(sizeof(sysconf_gni_peer_info_t) *
				       nranks);
	if (sysconf_peer_gni_info == NULL) {
		fprintf(stderr,"SYSCONF(%s,%d): malloc of %d bytes failed\n",
			sysconf_my_gni_info.cname,
			sysconf_my_gni_info.rank,
			sizeof(sysconf_gni_peer_info_t) * nranks);
		MPI_Abort(MPI_COMM_WORLD,-1);
	}

	for (i=0; i<nranks; i++) {
		assert(tmp_info[i].rank < nranks);
		sysconf_peer_gni_info[tmp_info[i].rank].gni_addr =
			tmp_info[i].gni_addr;
		strcpy(sysconf_peer_gni_info[tmp_info[i].rank].cname,
		       tmp_info[i].cname);
	}

	free(tmp_info);
	ugni_initialized = 1;

	return 0;
}

int ugni_finalize(void)
{
	gni_return_t status;

	/*
	 * wasn't initialized, nothing to do
	 */

	fprintf(stderr,"CALLING UGNI_FINALIZE\n");
	if (!ugni_initialized)
		return -1;

	status = GNI_CqDestroy(gni_nic.tx_cq);
	if (status != GNI_RC_SUCCESS) {
	}

	status = GNI_CqDestroy(gni_nic.rx_cq);
	if (status != GNI_RC_SUCCESS) {
	}

	status = GNI_CdmDestroy(gni_nic.cdm);
	if (status != GNI_RC_SUCCESS) {
	}

#if 0
	free(sysconf_peer_gni_info);
	sysconf_peer_gni_info = NULL;
#endif

	return 0;
}

static int _get_cname(char *cname, int len)
{
	int ret = 0;
	FILE *fd=NULL;
	char filename[1024];

	fd = fopen("/proc/cray_xt/cname","r");
	if (fd != NULL) {
		fscanf(fd, "%s", filename);
		if (strlen(filename) < len) {
			strcpy(cname, filename);
		} else {
			ret = -1;
		}
		fclose(fd);
	} else {
		ret = -1;
	}

	return ret;
}

static int _get_ptag(uint8_t *out_ptag)
{
    /* TODO no need for tmp */
    char *ptr;
    uint8_t tmp_ptag;

	if (NULL == (ptr = getenv("PMI_GNI_PTAG"))) {
		return -1;
	}

	errno = 0;
	tmp_ptag = (uint8_t)strtoul (ptr, (char **)NULL, 10);
	if (0 != errno) {
		return -1;
	}

	*out_ptag = tmp_ptag;
	return 0;
}

static int _get_cookie (uint32_t *out_cookie)
{
    char *ptr;
    uint32_t tmp_cookie;

    if (NULL == (ptr = getenv("PMI_GNI_COOKIE"))) {
        return -1;
    }

    errno = 0;
    tmp_cookie = (uint32_t) strtoul (ptr, NULL, 10);
    if (0 != errno) {
        /* TODO add err msg - better rc? */
        return -1;
    }

    *out_cookie = tmp_cookie;

    return 0;
}
