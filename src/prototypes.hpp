extern void Init_variables();
extern double lagrange(double*,double,int,int);
extern double lagrange_prime(double*,double,int,int);
extern void lagrange_matrix(double*,double*,double*,int,int);
extern void lagrange_prime_matrix(double*,double*,double*,int,int);
extern void gauss_legendre(double, double, int, double*, double*);
extern void flux_points(double*, double*, int);
extern void solution_points(double *, int);
extern void ader_matrix(double*, double*, double*, int);
extern void integral_matrix(double*, double*, double*, double*, double*, int);
extern void inverse(double*, double*, int);

extern void Build_mesh();
extern void build_faces(double*, double*, int, double);
extern void build_centers(double*, double*, int);


extern void prim_to_cons(double*, double*, int);
extern void cons_to_prim(double*, double*, int);
extern void fluxes(double*, double*, double*, int, int, int, int);
extern double sound_speed(double, double);
extern void riemann_llf(double*, double*, double*, int, int, int, int);
extern void riemann_hllc(double*, double*, double*, int, int, int, int);
extern void riemann_solver_x();
extern void riemann_solver_y();
extern void riemann_solver_z();

extern void transform_cv_to_sp(double*, double*);
extern void transform_sp_to_cv(double*, double*);
extern void transform_sp_to_fp_x(double*, double*, int);
extern void transform_sp_to_fp_y(double*, double*, int);
extern void transform_sp_to_fp_z(double*, double*, int);
extern void derive_fp_x_to_sp(double*, double*);
extern void derive_fp_y_to_sp(double*, double*);
extern void derive_fp_z_to_sp(double*, double*);
extern void face_integral_x(double*, double*, int);
extern void face_integral_y(double*, double*, int);
extern void face_integral_z(double*, double*, int);
extern void volume_integral(double*, double*);

extern void ader_subupdate(double*, double*, double*);
extern void ader_update(double*, double*);
extern void fv_update(double *, double *, int);

extern void Initial_Conditions( );
extern void Split_Domain( );
extern void Boundary_Conditions(double*);
extern void Store_boundaries(double*);
extern void store_boundaries_x(double*);
extern void store_boundaries_y(double*);
extern void store_boundaries_z(double*);
extern void apply_boundaries_x(double*);
extern void apply_boundaries_y(double*);
extern void apply_boundaries_z(double*);

//Time step
extern void compute_dt();

//Output functions
extern void Write_fields(int);
extern void Write_array(double*,int,char*);
extern void Write(int);
extern void Write_edges();

//Comms
extern void Build_comms();
extern void Exec_comms(double*);
extern void fill_comms_x(double*, double*, int);
extern void fill_comms_y(double*, double*, int);
extern void fill_comms_z(double*, double*, int);

extern void detect_troubles();
extern void godunov_2O();
extern void fv_godunov2O_update(double*, double *, int);
extern void fv_corrected_update(double*, double *, int);

extern void extrema(double*, double*, double*, int, double);
extern void smooth_extrema_x(double*, double*, int);
extern void smooth_extrema_y(double*, double*, int);
extern void smooth_extrema_z(double*, double*, int);

extern void finish();