/*
 * header file for ugni util related
 * struct definitions and function prototypes
 * to use uGNI direct in sysconfidence
 */

#include <stdio.h>
#include <stdint.h>
#include <fcntl.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <strings.h>
#include <assert.h>
#include <malloc.h>
#include <sched.h>
#include "gni_pub.h"
#include "pmi.h"

/*
 * typedefs and macros
 */

/*
 *note about cnames - they are of the format
 *
 * cA-BcCsXnY
 *
 * Where A,B indicate the col/row of the cabinet
 * in the system - this will at most need 2 characters
 * each,
 *
 * C is the chasis within a cabinet - this will range
 * from 0-2
 *
 * X is the slot in a chasis - this ranges from 0-15
 * Y is the node on a board - this ranges from 0-3
 *
 * so bottom line is 32 chars is enough for the cname
 */

typedef struct sysconf_gni_peer_info {
	int rank;
	uint32_t gni_addr;
	char cname[32];
} sysconf_gni_peer_info_t;


typedef struct sysconf_gni_nic {
	gni_cdm_handle_t cdm;
	gni_nic_handle_t nic;
	gni_cq_handle_t tx_cq;
	gni_cq_handle_t rx_cq;
} sysconf_gni_nic_t;

/*
 * globals
 */

extern sysconf_gni_nic_t gni_nic;
extern sysconf_gni_peer_info_t *sysconf_peer_gni_info;
extern sysconf_gni_peer_info_t sysconf_my_gni_info;
extern uint16_t sysconf_ugni_dlvr_mode;

/*
 * prototypes
 */

/**
 * @brief initialize ugni direct communication layer
 *
 * @return 0 on success, -1 on failure (actually aborts)
 */
int ugni_init(void);

/**
 * @brief finalize ugni direct communication layer and free up resources
 *        associated with the ugni comm package
 *
 * @return 0 on success, -1 on failure
 */
int ugni_finalize(void);

