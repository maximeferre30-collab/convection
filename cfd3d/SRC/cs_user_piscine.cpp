/*============================================================================
 * Niveau 3 : piscine de combustible usé en 3D – Code_Saturne 9.0
 *
 * Mêmes hypothèses que le 0D/2D (paquet Python piscine/) :
 *  - piscine 12 x 10 x 12,5 m, maillage cartésien interne (pas de Salome) ;
 *  - Boussinesq (idilat = 0) avec β(T) sécant : la flottabilité reproduit
 *    exactement ρ(T_init) - ρ(T) du fit utilisé en 0D ;
 *  - turbulence k-ω SST avec production par flottabilité ;
 *  - racks = bloc poreux : frein de Darcy implicite + puissance de Way-Wigner ;
 *  - surface libre = symétrie (glissement) + flux de pertes flux_pertes(T) ;
 *  - fond et parois : paroi lisse adiabatique.
 *
 * Variables d'environnement (lues par lancer.sh) :
 *   PISCINE_REPARTITION  uniforme | bord | coin   (assemblages chauds x3)
 *   PISCINE_NX           cellules selon x (y et z en proportion), défaut 48
 *   PISCINE_DUREE_H      durée physique simulée (h), défaut 6
 *   PISCINE_SORTIE_H     période des instantanés CSV (h), défaut 1
 *============================================================================*/

#include "base/cs_defs.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "cs_headers.h"

BEGIN_C_DECLS

/*----------------------------------------------------------------------------
 * Paramètres physiques (identiques à piscine/modele0d.py, cas de référence)
 *----------------------------------------------------------------------------*/

static const double LX = 12.0, LY = 10.0, H = 12.5;          /* m */
static const double RX0 = 1.0, RX1 = 11.0;                    /* racks (m) */
static const double RY0 = 1.0, RY1 = 9.0;
static const double RZ0 = 0.3, RZ1 = 4.5;
static const double T_INIT = 30.0, T_AIR = 30.0, HR = 0.5;    /* °C */
static const double P0 = 1000e6;                               /* W */
static const double T_IRR = 3*365.25*86400., T_REFROID = 30*86400.;
static const double SIGMA_RACK = 0.2;                          /* 1/s */
static const double CP = 4186.0, LAMBDA = 0.62;                /* J/kg/K, W/m/K */
static const double TK = 273.15;

static double
_env(const char *nom, double defaut)
{
  const char *s = getenv(nom);
  return (s != nullptr && s[0] != '\0') ? atof(s) : defaut;
}

static const char *
_repartition(void)
{
  const char *s = getenv("PISCINE_REPARTITION");
  return (s != nullptr && s[0] != '\0') ? s : "uniforme";
}

static double _rho(double T) { return 1002.305 - 0.13347*T - 3.0766e-3*T*T; }

static double _psat(double T)
{
  return 133.322*pow(10., 8.07131 - 1730.63/(233.426 + T));
}

/* Flux perdu par la surface (W/m2) : évaporation ASHRAE, convection, rayonnement */
static double
_flux_pertes(double T)
{
  double dp = fmax(_psat(T) - HR*_psat(T_AIR), 0.)/1000.;
  double q_evap = dp*0.089*1000.;
  double dT = fmax(T - T_AIR, 0.);
  double q_conv = 1.52*pow(dT, 4./3.);
  double q_rad = 0.95*5.670e-8*(pow(T + TK, 4) - pow(T_AIR + TK, 4));
  return q_evap + q_conv + q_rad;
}

static double
_puissance(double t)
{
  double ts = T_REFROID + t;
  return 0.0622*P0*(pow(ts, -0.2) - pow(ts + T_IRR, -0.2));
}

static bool
_dans_racks(const cs_real_t x[3])
{
  return    x[0] >= RX0 && x[0] <= RX1 && x[1] >= RY0 && x[1] <= RY1
         && x[2] >= RZ0 && x[2] <= RZ1;
}

/* Poids de puissance : tiers chaud x3 selon x ("bord", comme le 2D)
   ou coin x/y ("coin", asymétrie 3D) */
static double
_poids(const cs_real_t x[3])
{
  if (!_dans_racks(x))
    return 0.;
  const char *rep = _repartition();
  bool x_chaud = x[0] < RX0 + (RX1 - RX0)/3.;
  bool y_chaud = x[1] < RY0 + (RY1 - RY0)/3.;
  if (strcmp(rep, "bord") == 0 && x_chaud)
    return 3.;
  if (strcmp(rep, "coin") == 0 && x_chaud && y_chaud)
    return 3.;
  return 1.;
}

/*----------------------------------------------------------------------------
 * Maillage cartésien
 *----------------------------------------------------------------------------*/

void
cs_user_mesh_cartesian_define(void)
{
  int nx = (int)_env("PISCINE_NX", 48);
  int nxyz[3] = {nx, (int)lround(nx*LY/LX), (int)lround(nx*H/LX)};
  cs_real_t xyz[6] = {0., 0., 0., LX, LY, H};
  cs_mesh_cartesian_define_simple(nullptr, nxyz, xyz);
}

/*----------------------------------------------------------------------------
 * Zones de bord
 *----------------------------------------------------------------------------*/

void
cs_user_zones(void)
{
  cs_boundary_zone_define("surface", "Z1", 0);
  cs_boundary_zone_define("parois", "X0 or X1 or Y0 or Y1 or Z0", 0);
}

/*----------------------------------------------------------------------------
 * Choix des modèles
 *----------------------------------------------------------------------------*/

void
cs_user_model(void)
{
  cs_thermal_model_t *th = cs_get_glob_thermal_model();
  th->thermal_variable = CS_THERMAL_MODEL_TEMPERATURE;
  th->temperature_scale = CS_TEMPERATURE_SCALE_KELVIN;

  cs_turb_model_t *turb = cs_get_glob_turb_model();
  turb->model = CS_TURB_K_OMEGA;   /* k-omega SST */

  cs_time_step_options_t *to = cs_get_glob_time_step_options();
  to->idtvar = CS_TIME_STEP_ADAPTIVE;
}

/*----------------------------------------------------------------------------
 * Paramètres numériques et propriétés de référence
 *----------------------------------------------------------------------------*/

void
cs_user_parameters(cs_domain_t *domain)
{
  cs_get_glob_velocity_pressure_model()->idilat = 0;   /* Boussinesq */

  cs_physical_constants_t *pc = cs_get_glob_physical_constants();
  pc->gravity[0] = 0.; pc->gravity[1] = 0.; pc->gravity[2] = -9.81;

  cs_fluid_properties_t *fp = cs_get_glob_fluid_properties();
  fp->ro0 = _rho(T_INIT);
  fp->viscl0 = 7.97e-4;            /* recalculée par cellule (Vogel) */
  fp->cp0 = CP;
  fp->t0 = T_INIT + TK;            /* référence de Boussinesq */
  fp->p0 = 101325.;
  fp->irovar = 0;
  fp->ivivar = 1;
  fp->icp = -1;

  cs_field_set_key_double(CS_F_(t), cs_field_key_id("diffusivity_ref"), LAMBDA);

  cs_get_glob_turb_rans_model()->has_buoyant_term = 1;

  cs_time_step_options_t *to = cs_get_glob_time_step_options();
  /* Courant max 4 : schéma implicite, et le transitoire étudié (heures) est
     bien plus lent que le temps de transit d'une maille (~2 s). */
  to->coumax = 4.;
  to->dtmax = 10.;
  domain->time_step->dt_ref = 0.5;
  domain->time_step->t_max = 3600.*_env("PISCINE_DUREE_H", 6.);

  cs_log_iteration_set_interval(50);
}

/*----------------------------------------------------------------------------
 * Sorties EnSight pour ParaView (toutes les 30 min de temps physique)
 *----------------------------------------------------------------------------*/

void
cs_user_postprocess_writers(void)
{
  cs_post_define_writer(CS_POST_WRITER_DEFAULT, "results", "postprocessing",
                        "EnSight Gold", "", FVM_WRITER_FIXED_MESH,
                        false, true, -1, 1800.);
}

/*----------------------------------------------------------------------------
 * Initialisation : eau au repos à T_INIT
 *----------------------------------------------------------------------------*/

void
cs_user_initialization(cs_domain_t *domain)
{
  if (domain->time_step->nt_prev > 0)
    return;  /* reprise */
  cs_real_t *T = CS_F_(t)->val;
  for (cs_lnum_t c = 0; c < domain->mesh->n_cells; c++)
    T[c] = T_INIT + TK;
}

/*----------------------------------------------------------------------------
 * Propriétés variables : μ(T) (Vogel) et β(T) sécant pour Boussinesq
 *----------------------------------------------------------------------------*/

void
cs_user_physical_properties(cs_domain_t *domain)
{
  const cs_real_t *T = CS_F_(t)->val;
  cs_real_t *mu = CS_F_(mu)->val;
  cs_real_t *beta = cs_field_by_name("thermal_expansion")->val;
  const double rho0 = _rho(T_INIT);

  for (cs_lnum_t c = 0; c < domain->mesh->n_cells; c++) {
    double Tc = T[c] - TK;
    mu[c] = 2.414e-5*pow(10., 247.8/(T[c] - 140.));
    double dT = Tc - T_INIT;
    /* ρ0 β (T - T0) = ρ(T0) - ρ(T) exactement (dérivée si dT -> 0) */
    beta[c] = (fabs(dT) > 1e-6) ? (rho0 - _rho(Tc))/(rho0*dT)
                                : (0.13347 + 2*3.0766e-3*Tc)/rho0;
  }
}

/*----------------------------------------------------------------------------
 * Termes sources : frein des racks (vitesse), puissance résiduelle (T)
 *----------------------------------------------------------------------------*/

void
cs_user_source_terms(cs_domain_t  *domain,
                     int           f_id,
                     cs_real_t    *st_exp,
                     cs_real_t    *st_imp)
{
  const cs_field_t *f = cs_field_by_id(f_id);
  const cs_lnum_t n_cells = domain->mesh->n_cells;
  const cs_real_t *vol = domain->mesh_quantities->cell_vol;
  const cs_real_3_t *xc = domain->mesh_quantities->cell_cen;
  const cs_real_t *rho = CS_F_(rho)->val;

  if (f == CS_F_(vel)) {
    cs_real_33_t *imp = (cs_real_33_t *)st_imp;
    for (cs_lnum_t c = 0; c < n_cells; c++) {
      if (_dans_racks(xc[c]))
        for (int k = 0; k < 3; k++)
          imp[c][k][k] -= vol[c]*rho[c]*SIGMA_RACK;
    }
  }
  else if (f == CS_F_(t)) {
    double somme = 0.;
    for (cs_lnum_t c = 0; c < n_cells; c++)
      somme += _poids(xc[c])*vol[c];
    cs_parall_sum(1, CS_DOUBLE, &somme);
    const double P = _puissance(domain->time_step->t_cur);
    for (cs_lnum_t c = 0; c < n_cells; c++)
      st_exp[c] += P*_poids(xc[c])*vol[c]/somme;   /* W (Cp hors dérivée) */
  }
}

/*----------------------------------------------------------------------------
 * Conditions aux limites
 *----------------------------------------------------------------------------*/

void
cs_user_boundary_conditions(cs_domain_t  *domain,
                            int           bc_type[])
{
  const cs_lnum_t *b_face_cells = domain->mesh->b_face_cells;
  const cs_real_t *T = CS_F_(t)->val;
  cs_field_t *th = CS_F_(t);

  const cs_zone_t *zp = cs_boundary_zone_by_name("parois");
  for (cs_lnum_t i = 0; i < zp->n_elts; i++)
    bc_type[zp->elt_ids[i]] = CS_SMOOTHWALL;   /* flux thermique nul par défaut */

  const cs_zone_t *zs = cs_boundary_zone_by_name("surface");
  for (cs_lnum_t i = 0; i < zs->n_elts; i++) {
    cs_lnum_t f_id = zs->elt_ids[i];
    bc_type[f_id] = CS_SYMMETRY;
    th->bc_coeffs->icodcl[f_id] = 3;
    /* flux sortant positif (W/m2), évalué sur la maille adjacente */
    th->bc_coeffs->rcodcl3[f_id] = _flux_pertes(T[b_face_cells[f_id]] - TK);
  }
}

/*----------------------------------------------------------------------------
 * Séries temporelles (bilan d'énergie, points chauds) et instantanés CSV
 *----------------------------------------------------------------------------*/

void
cs_user_extra_operations(cs_domain_t *domain)
{
  const cs_mesh_t *m = domain->mesh;
  const cs_mesh_quantities_t *mq = domain->mesh_quantities;
  const cs_time_step_t *ts = domain->time_step;
  const cs_real_t *T = CS_F_(t)->val;
  const cs_real_3_t *u = (const cs_real_3_t *)CS_F_(vel)->val;
  const cs_real_t *vol = mq->cell_vol;
  const cs_real_3_t *xc = mq->cell_cen;

  static double t_serie = -1e30, t_snap = 0.;
  static bool entete = false;
  const double periode_snap = 3600.*_env("PISCINE_SORTIE_H", 1.);
  bool fin = (ts->t_cur >= ts->t_max - 1e-9) || (ts->nt_cur >= ts->nt_max);

  if (ts->t_cur >= t_serie + 300. || fin) {   /* toutes les 5 min */
    t_serie = ts->t_cur;
    double s[4] = {0., 0., 0., 0.};   /* ∫T dV, V, ∫T_surf dA, A */
    double mx[2] = {-1e30, 0.};       /* T_max, |u|_max */
    for (cs_lnum_t c = 0; c < m->n_cells; c++) {
      s[0] += T[c]*vol[c];
      s[1] += vol[c];
      mx[0] = fmax(mx[0], T[c]);
      mx[1] = fmax(mx[1], cs_math_3_norm(u[c]));
    }
    double q_surf = 0.;
    const cs_zone_t *zs = cs_boundary_zone_by_name("surface");
    for (cs_lnum_t i = 0; i < zs->n_elts; i++) {
      cs_lnum_t f_id = zs->elt_ids[i];
      double Tf = T[m->b_face_cells[f_id]];
      double a = mq->b_face_surf[f_id];
      s[2] += Tf*a;
      s[3] += a;
      q_surf += _flux_pertes(Tf - TK)*a;
    }
    cs_parall_sum(4, CS_DOUBLE, s);
    cs_parall_sum(1, CS_DOUBLE, &q_surf);
    cs_parall_max(2, CS_DOUBLE, mx);

    if (cs_glob_rank_id <= 0) {
      FILE *fp = fopen("piscine_series.csv", entete ? "a" : "w");
      if (!entete)
        fprintf(fp, "t_s,T_moy_C,T_max_C,T_surf_C,v_max_m_s,P_W,pertes_W,"
                    "energie_J\n");
      entete = true;
      double rc = _rho(T_INIT)*CP;
      fprintf(fp, "%.1f,%.5f,%.5f,%.5f,%.5e,%.6e,%.6e,%.10e\n",
              ts->t_cur, s[0]/s[1] - TK, mx[0] - TK, s[2]/s[3] - TK, mx[1],
              _puissance(ts->t_cur), q_surf, rc*s[0]);
      fclose(fp);
    }
  }

  /* Instantanés : un fichier par rang (x, y, z, T, u, v, w), fusionnés en Python */
  if (ts->t_cur >= t_snap - 1e-9 || fin) {
    char nom[64];
    snprintf(nom, 64, "champ_t%06.0fs_r%03d.csv", ts->t_cur,
             cs_glob_rank_id < 0 ? 0 : cs_glob_rank_id);
    FILE *fp = fopen(nom, "w");
    fprintf(fp, "x,y,z,T,u,v,w\n");
    for (cs_lnum_t c = 0; c < m->n_cells; c++)
      fprintf(fp, "%.4f,%.4f,%.4f,%.4f,%.4e,%.4e,%.4e\n",
              xc[c][0], xc[c][1], xc[c][2], T[c] - TK, u[c][0], u[c][1], u[c][2]);
    fclose(fp);
    while (t_snap <= ts->t_cur + 1e-9)
      t_snap += periode_snap;
  }
}

END_C_DECLS
