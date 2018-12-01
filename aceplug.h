#ifndef ___PLUGIN
#define ___PLUGIN 1

#define ACEPLUG_VERSION 306

#ifndef __VECTOR_TYPES_H__

typedef struct {
  float x; float y; float z; float w;
} float4;

typedef struct {
  double x; double y; double z; double w;
} double4;

#endif


  typedef struct aceplug_pdb_t {
		double4 pos;
    /*
     *                        definition                         columns
     */
    char record[8];  /**< Record name: "ATOM  " or "HETATM"   (col:  1 -  6) */
    char serial[8];  /**< Atom serial number.                 (col:  7 - 11) */
    char name[8];    /**< Atom name.                          (col: 13 - 16) */
    char altLoc[4];  /**< Alternate location identifier.      (col:      17) */
    char resName[4]; /**< Residue name.                       (col: 18 - 20) */
    char chainID[4]; /**< Chain identifier.                   (col:      22) */
    char resSeq[8];  /**< Residue sequence number.            (col: 23 - 26) */
    char iCode[4];   /**< Code for insertion of residues.     (col:      27) */
    float occupancy; /**< Occupancy.                          (col: 55 - 60) */
    float tempFactor;/**< Temperature factor.                 (col: 61 - 66) */
    char segID[8];   /**< Segment identifier, left-justified. (col: 73 - 76) */
    char element[4]; /**< Element symbol, right-justified.    (col: 77 - 78) */
    char charge[4];  /**< Charge on the atom.                 (col: 79 - 80) */
  } aceplug_pdb_t;


  /* Time unit in femtoseconds */
#define   TIMEFACTOR   48.88821 
  /* Coulomb's constant for electrostatics, units kcal*A/mol/e^2 */
#define MD_COULOMB  332.0636
  /* Boltzman's constant for temperature, units kcal/A/K */
#define MD_BOLTZMAN  0.001987191

#define MD_KCAL_TO_KJ  4.1868


#define XCONCAT(x,y,z) x##y##z
#define XCONCAT2(x,y,z) XCONCAT(x,y,z)
#define _FUNCTION(prefix,suffix) XCONCAT2(prefix,_,suffix)

typedef enum {
	ACEPLUG_OK=0,
	ACEPLUG_INVALID,
	ACEPLUG_ERR,
	ACEPLUG_UNIMPL
} aceplug_err_t;

typedef enum {
	ENERGY_PE,
	ENERGY_KE,
	TEMPERATURE,
	ENERGY_ELEC,
	ENERGY_VDW,
	ENERGY_EXTERNAL,
	ENERGY_MAX=32
	
} aceplug_energy_t;





typedef struct aceplug_sim_t {
	int version;
	double4 *pos;
	float4 *frc;
	float4 *vel;
	float  *mass;
	float  *charge;
	float4 *vdw; // epsion, sigma, epsilon14 sigma14 
	int    natoms;
	unsigned long step;
	void *privdata;
	float4 box;
	float timestep_fs;

	int ensemble_rank;
	int ensemble_size;

// API call function pointers
	aceplug_err_t (*plugin_get_temperature)( double *temp );
	aceplug_err_t (*plugin_set_temperature)( double val );
	aceplug_err_t (*plugin_register_tcl_function)( char *command, /*Tcl_ObjCmdProc*/void *func, void *data ) ;

	aceplug_err_t (*plugin_load_positions)( void );
	aceplug_err_t (*plugin_update_positions)( void );

	aceplug_err_t (*plugin_load_velocities)( void );
	aceplug_err_t (*plugin_update_velocities)( void );

	aceplug_err_t (*plugin_load_forces)( void );
	aceplug_err_t (*plugin_update_forces)( void );


	aceplug_err_t (*plugin_set_pairlist) ( int *indices, int N );
	aceplug_err_t (*plugin_get_pairlist) ( int **pairlist, size_t *offset, int *N );
	aceplug_err_t (*plugin_add_energy) ( aceplug_energy_t term, double energy );
	aceplug_err_t (*plugin_add_virial) (  double virial[3][3] );

	aceplug_err_t (*plugin_ensemble_send) ( int peer, void *data, size_t len );
	aceplug_err_t (*plugin_ensemble_recv) ( int peer, void *data, size_t len );

	aceplug_err_t (*plugin_get_energy_temp)( double *val );
	aceplug_err_t (*plugin_vdw_param)( int iatom, int jatom, int s14, double *A, double *B );

	aceplug_err_t (*plugin_read_pdb)( char *filename, struct aceplug_pdb_t **pdb, int *N );
	aceplug_err_t (*plugin_write_pdb)( char *filename, struct aceplug_pdb_t *pdb, int N );

	int (*plugin_version)( void);

	void (*plugin_scale_velocities)( double factor );

	aceplug_err_t (*plugin_exchange_system)( int peer );
	aceplug_err_t (*plugin_load_vdwtable)( void **table, size_t *len );
	aceplug_err_t (*plugin_update_vdwtable)( void *table );
	aceplug_err_t (*plugin_load_charges)( void **table, size_t *len );
	aceplug_err_t (*plugin_update_charges)( void *table );
	aceplug_err_t (*plugin_reevaluate_energies)( double *energies );

	aceplug_err_t (*plugin_ensemble_alltoall) ( void *sendbuf, size_t len, void *recvbuf );
	aceplug_err_t (*plugin_ensemble_allgather) ( void *sendbuf, size_t len, void *recvbuf );
} aceplug_sim_t ;

#ifdef __cplusplus
extern "C" {
#endif

// API Calls

aceplug_err_t aceplug_init( aceplug_sim_t *s, int argc, char **argkey, char **argval );
aceplug_err_t aceplug_calcforces( aceplug_sim_t *s );
aceplug_err_t aceplug_endstep( aceplug_sim_t *s );
aceplug_err_t aceplug_terminate( aceplug_sim_t *s );

#ifdef __cplusplus
}
#endif


#endif



