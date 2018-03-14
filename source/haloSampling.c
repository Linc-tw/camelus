

  /*******************************************************
   **  haloSampling.c					**
   **  Version 2018.03.13				**
   **							**
   **  Copyright (C) 2018 - Chieh-An Lin		**
   **  GNU GPLv3 - https://www.gnu.org/licenses/	**
   *******************************************************/


#include "haloSampling.h"


//----------------------------------------------------------------------
//-- Functions related to halo_t, halo_node, halo_list

void set_halo_t(cosmo_hm *chPar, peak_param *pkPar, halo_t *h, double pos[2], double z, double ww, double M, error **err)
{
  testErrorRet(h==NULL, peak_null, "halo_t *h is NULL", *err, __LINE__,);
  
  if (pos != NULL) {
    h->pos[0] = pos[0];
    h->pos[1] = pos[1];
    
    if (pkPar->field > 0) {
      h->sinCosDEC[0] = sin(h->pos[1]);
      h->sinCosDEC[1] = cos(h->pos[1]);
    }
  }
  
  h->z      = z;
  double a  = 1.0 / (1.0 + z);
  h->a      = a;
  if (ww < 0.0) {
    h->w    = w(chPar->cosmo, a, 0, err);
    forwardError(*err, __LINE__,); //-- wOmegar = 0
  }
  else h->w = ww;
  h->M      = M;
  double c  = concentration(chPar, M, a, err);  forwardError(*err, __LINE__,);
  h->c      = c;
  h->c_sq   = SQ(c);
  h->f      = 1.0 / (log(1.0 + c) - c/(1.0 + c));
  
  //-- r_vir is halo's physical size
  //-- M = (4 pi / 3) * rho_vir(z) * r_vir^3
  //--   = (4 pi / 3) * (rho_vir(z) / rho_m(z)) * rho_m(z) * r_vir^3
  //--   = (4 pi / 3) * Delta_vir(z) * rho_m(0) * (1+z)^3 * r_vir^3
  //--   = (4 pi / 3) * Delta_vir(z) * rho_crit(0) * Omega_m * (r_vir / a)^3
  h->r_vir        = a * cbrt(M / (FOUR_PI_OVER_THREE * Delta_vir(chPar, a) * CRITICAL_DENSITY * chPar->cosmo->Omega_m));
  h->r_vir_sq     = SQ(h->r_vir);
  //-- theta_vir = r_vir / D_A
  //-- given by angular diameter distance D_A = f_K / (1+z)
  double D_l      = a * f_K(chPar->cosmo, h->w, err); forwardError(*err, __LINE__,);
  h->theta_vir    = h->r_vir / D_l;
  h->theta_vir_sq = SQ(h->theta_vir);
  
  //-- factor = FOUR_PI_G_OVER_C2 * (D_l M f_NFW c_NFW^2 / 2 pi r_vir^2)
  //-- See rayTracing.c for details
  h->factor = FOUR_PI_G_OVER_C2 * D_l * M * h->f * h->c_sq / (TWO_PI * h->r_vir_sq);
  return;
}

halo_node *initialize_halo_node(error **err)
{
  halo_node *hNode = (halo_node*)malloc_err(sizeof(halo_node), err); forwardError(*err, __LINE__, NULL);
  hNode->h         = (halo_t*)malloc_err(sizeof(halo_t), err);       forwardError(*err, __LINE__, NULL);
  hNode->next      = NULL;
  return hNode;
}

halo_list *initialize_halo_list(error **err)
{
  halo_list *hList = (halo_list*)malloc_err(sizeof(halo_list), err); forwardError(*err, __LINE__, NULL);
  hList->length    = 0;
  hList->size      = 0;
  hList->first     = NULL;
  hList->last      = NULL;
  return hList;
}

void free_halo_list(halo_list *hList)
{
  halo_node *hNode;
  if (hList) {
    while (hList->first != NULL) {
      hNode        = hList->first;
      hList->first = hNode->next;
      if (hNode->h) {free(hNode->h); hNode->h = NULL;}
      free(hNode); hNode = NULL;
    }
    free(hList); hList = NULL;
  }
  return;
}

void reset_halo_list(halo_list *hList)
{
  hList->size = 0;
  return;
}

void append_halo_list(cosmo_hm* chPar, peak_param *pkPar, halo_list* hList, double pos[2], double z, double ww, double M, error** err)
{
  halo_t *h;
  if (hList->length == 0) {
    hList->first = initialize_halo_node(err);      forwardError(*err, __LINE__,);
    hList->last  = hList->first;
    hList->length++;
  }
  else if (hList->length == hList->size) {
    hList->last->next = initialize_halo_node(err); forwardError(*err, __LINE__,);
    hList->last       = hList->last->next;
    hList->length++;
  }
  else if (hList->size == 0) hList->last = hList->first;
  else                       hList->last = hList->last->next;
  
  h = hList->last->h;
  set_halo_t(chPar, pkPar, h, pos, z, ww, M, err); forwardError(*err, __LINE__,);
  hList->size++;
  return;
}

//----------------------------------------------------------------------
//-- Functions related to halo_map

halo_map *initialize_halo_map(int N1, int N2, double theta_pix, error **err)
{
  halo_map *hMap      = (halo_map*)malloc_err(sizeof(halo_map), err);                    forwardError(*err, __LINE__, NULL);
  hMap->N1            = N1;
  hMap->N2            = N2;
  hMap->length        = N1 * N2;
  hMap->total         = 0;
  hMap->theta_pix     = theta_pix;
  hMap->theta_pix_inv = 1.0 / theta_pix;
  
  hMap->map           = (halo_list**)malloc_err(hMap->length * sizeof(halo_list*), err); forwardError(*err, __LINE__, NULL);
  int i;
  for (i=0; i<hMap->length; i++) hMap->map[i] = initialize_halo_list(err);
  return hMap;
}

void free_halo_map(halo_map *hMap)
{
  int i;
  if (hMap) {
    if (hMap->map) {
      for (i=0; i<hMap->length; i++) {free_halo_list(hMap->map[i]); hMap->map[i] = NULL;}
      free(hMap->map); hMap->map = NULL;
    }
    free(hMap); hMap = NULL;
  }
  return;
}

void reset_halo_map(halo_map *hMap)
{
  hMap->total = 0;
  int i;
  for (i=0; i<hMap->length; i++) reset_halo_list(hMap->map[i]);
  return;
}

int append_halo_map(cosmo_hm *chPar, peak_param *pkPar, halo_map *hMap, double pos[2], double z, double ww, double M, error **err)
{
  double pos2[2];
  if (pkPar->field == 0) {
    pos2[0] = pos[0];
    pos2[1] = pos[1];
  }
  else {
    project(pos, pkPar->center, pos2);
    rotate(pos2, pkPar->rotAng, pos2);
    pos2[0] += pkPar->offset;
    pos2[1] += pkPar->offset;
  }
  
  int i = (int)(pos2[0] * hMap->theta_pix_inv);
  int j = (int)(pos2[1] * hMap->theta_pix_inv);
  if ((i < 0) || (j < 0) || (i >= hMap->N1) || (j >= hMap->N2)) return 0;
  append_halo_list(chPar, pkPar, hMap->map[i+j*hMap->N1], pos, z, ww, M, err); forwardError(*err, __LINE__, -1);
  hMap->total++;
  return 1;
}

void read_halo_map(char name[], cosmo_hm *chPar, peak_param *pkPar, halo_map *hMap, error **err)
{
  reset_halo_map(hMap);
  
  //-- Open
  FILE *file = fopen_err(name, "r", err); forwardError(*err, __LINE__,);
  if (pkPar->verbose < 2) {
    printf("Reading...\r"); fflush(stdout);
  }
  
  //-- Column indices, inHaloCatCol successively read as pos1, pos2, z, w, M
  int pos1Ind = pkPar->inHaloCatCol[0];
  int pos2Ind = pkPar->inHaloCatCol[1];
  int zInd    = pkPar->inHaloCatCol[2];
  int wInd    = pkPar->inHaloCatCol[3];
  int MInd    = pkPar->inHaloCatCol[4];
  
  int count  = 0;
  int count2 = 0;
  double factor = (pkPar->field == 0) ? ARCMIN_TO_RADIAN : DEGREE_TO_RADIAN;
  
  char kv[STRING_LENGTH_MAX][256], buffer[STRING_LENGTH_MAX];
  double pos[2], z, ww, M;
  
  //-- Read
  int c = fgetc(file);
  while (c != EOF) {
    if (c == (int)'#') fgets(buffer, STRING_LENGTH_MAX, file);
    else {
      ungetc(c, file);
      fgets(buffer, STRING_LENGTH_MAX, file);
      getKeyAndValues(buffer, kv);
      count++;
      
      pos[0] = atof(kv[pos1Ind]) * factor; //-- [rad]
      pos[1] = atof(kv[pos2Ind]) * factor; //-- [rad]
      z      = atof(kv[zInd]);
      ww     = (wInd < 0) ? -1.0 : atof(kv[wInd]);
      M      = atof(kv[MInd]);
      count2 += append_halo_map(chPar, pkPar, hMap, pos, z, ww, M, err); forwardError(*err, __LINE__,);
    }
    c = fgetc(file);
  }
  fclose(file);
  
  //-- Check
  testErrorRet(count2!=hMap->total, peak_match, "Number not matched after reading files", *err, __LINE__,);
  if (pkPar->verbose < 3) {
    printf("Read \"%s\"   \n", name);
    printf("Found %d halos, generated %d\n", count, count2);
  }
  else if (pkPar->verbose == 3) {
    printf("Read %d halos, ", count2); fflush(stdout);
  }
  return;
}

//----------------------------------------------------------------------
//-- Functions related to mass function

double massFct(cosmo_hm_params *cANDp, double mass, error **err)
{
  //-- nicaea2.5 calculates with comoving volume, nicaea2.7 with physical.
  //-- Here I change it into a comoving number density.
  //-- The factor ln(10) corrects dn/dlnM to dn/dlogM.
  double dn = LN_10 * CB(cANDp->a) * dn_dlnM(mass, (void*)cANDp, err); forwardError(*err, __LINE__, -1.0);
  return dn;
}

void fillMassFct(cosmo_hm_params *cANDp, sampler_t *hSamp, error **err)
{
  double *x    = hSamp->x;
  double *pdf  = hSamp->pdf;
  int i;
  
  //-- Fill pdf
  for (i=0; i<hSamp->length; i++) {
    pdf[i] = massFct(cANDp, x[i], err);
    forwardError(*err, __LINE__,);
  }
  
  //-- Precompute cdf, total pdf, <x>
  int setTotalToOne = 0; //-- This is to keep the total of pdf.
  set_sampler_t(hSamp, setTotalToOne);
  //-- After the above operation:
  //--   pdf is normalized
  //--   cdf is normalized
  //--   totPdf is not normalized
  //--   x_mean is not normalized
  return;
}

void outAsciiMassFct(char name[], cosmo_hm *chPar, peak_param *pkPar, double z, error **err)
{
  cosmo_hm_params *cANDp = (cosmo_hm_params*)malloc_err(sizeof(cosmo_hm_params), err); forwardError(*err, __LINE__,);
  cANDp->model = chPar;
  cANDp->asymptotic = 0;
  cANDp->a = 1.0/(1.0+z);
  
  FILE *file = fopen_err(name, "w", err); forwardError(*err, __LINE__,);
  
  fprintf(file, "# Halo mass function, by comoving volume\n");
  fprintf(file, "#\n");
  fprintf(file, "# Key-in parameters\n");
  fprintf(file, "# z = %f\n", z);
  fprintf(file, "#\n");
  fprintf(file, "# Halo model parameters = %s\n", pkPar->hmParPath);
  fprintf(file, "# Mass function model = %s\n", smassfct_t(chPar->massfct));
  fprintf(file, "#\n");
  fprintf(file, "# Peak parameters = %s\n", pkPar->hmParPath);
  fprintf(file, "# M_min = %8.2e, M_max = %8.2e [M_sol/h], dlogM = %g\n", pkPar->M_min, pkPar->M_max, pkPar->dlogM);
  fprintf(file, "#\n");
  outAsciiCosmoParam(file, chPar, pkPar);
  fprintf(file, "#\n");
  fprintf(file, "#           M      dn/dlogM\n");
  fprintf(file, "#    [M_sol/h]  [(Mpc/h)^-3]\n");
  
  double M, dn;
  for (M=pkPar->M_min; M<=pkPar->M_max; M*=pow10(0.001)) {
    dn = massFct(cANDp, M, err); forwardError(*err, __LINE__,);
    //dn *= pow(1.0 + z, 3.0); //-- Decomment this to switch to physical volume
    fprintf(file, "  %12.6e  %12.6e\n", M, dn);
  }
  
  free(cANDp);
  fclose(file);
  if (pkPar->verbose < 3) printf("Outputed \"%s\"\n", name);
  return;
}

void outFitsMassFct(char name[], cosmo_hm *chPar, peak_param *pkPar, double z, error **err)
{
#ifdef __CAMELUS_USE_FITS__
  cosmo_hm_params *cANDp = (cosmo_hm_params*)malloc_err(sizeof(cosmo_hm_params), err); forwardError(*err, __LINE__,);
  cANDp->model = chPar;
  cANDp->asymptotic = 0;
  cANDp->a = 1.0/(1.0+z);
  
  FITS_t *fits = initializeTableWriter_FITS_t(name);
  
  addColumn(fits, "M",    TFLOAT, "M_sol/h");
  addColumn(fits, "NOFM", TFLOAT, "(Mpc/h)^-3");
  updateComment(fits, "TTYPE2", "Comoving volume");
  
  double M, dn;
  float fBuff;
  
  for (M=pkPar->M_min; M<=pkPar->M_max; M*=pow10(0.001)) {
    dn = massFct(cANDp, M, err); forwardError(*err, __LINE__,);
    //dn *= pow(1.0 + z, 3.0); //-- Decomment this to switch to physical volume
    fBuff = (float)M;  writeTableColumn(fits, 0, 1, &fBuff);
    fBuff = (float)dn; writeTableColumn(fits, 1, 1, &fBuff);
    nextRow(fits);
  }
  
  addKeyword(fits, TDOUBLE, "Z",       &z,                         "[-] Redshift");
  
  outFitsCosmoParam(fits, chPar, pkPar);
  
  addLineSpread(fits);
  addKeyword(fits, TSTRING, "HMPARAM", pkPar->hmParPath,           "[-] Path of the halo model parameter file");
  addKeyword(fits, TSTRING, "MASSFCT", smassfct_t(chPar->massfct), "[-] Mass function model");
  
  addLineSpread(fits);
  addKeyword(fits, TSTRING, "PKPARAM", pkPar->pkParPath,           "[-] Path of the peak parameter file");
  sprintf(name, "%8.2e", pkPar->M_min);
  addKeyword(fits, TSTRING, "MMIN",    &name,                      "[M_sol/h] Minimum halo mass");
  sprintf(name, "%8.2e", pkPar->M_max);
  addKeyword(fits, TSTRING, "MMAX",    &name,                      "[M_sol/h] Maximum halo mass");
  addKeyword(fits, TDOUBLE, "DLOGM",   &pkPar->dlogM,              "[-] Halo mass binwidth");
  
  free(cANDp);
  free_FITS_t(fits);
  if (pkPar->verbose < 3) printf("Outputed \"%s\"\n", name);
#endif
  return;
}

//----------------------------------------------------------------------
//-- Functions related to halo output

void outAscii_halo_map(FILE *file, peak_param *pkPar, halo_map *hMap)
{
  fprintf(file, "# Number of halos = %d\n", hMap->total);
  fprintf(file, "#\n");
  if (pkPar->field == 0) {
    fprintf(file, "#   theta_x    theta_y       z         w          M\n");
    fprintf(file, "#  [arcmin]   [arcmin]      [-]   [Mpc/h]  [M_sol/h]\n");
  }
  else {
    fprintf(file, "#       RA        DEC        z         w          M\n");
    fprintf(file, "#     [deg]      [deg]      [-]   [Mpc/h]  [M_sol/h]\n");
  }
  
  halo_list *hList;
  halo_node *hNode;
  halo_t *h;
  int i, j;
  if (pkPar->field == 0) {
    for (i=0; i<hMap->length; i++) {
      hList = hMap->map[i];
      for (j=0, hNode=hList->first; j<hList->size; j++, hNode=hNode->next) {
	h = hNode->h;
	fprintf(file, "  %9.3f  %9.3f  %7.5f  %8.3f  %9.3e\n", h->pos[0]*RADIAN_TO_ARCMIN, h->pos[1]*RADIAN_TO_ARCMIN, h->z, h->w, h->M);
      }
    }
  }
  else {
    for (i=0; i<hMap->length; i++) {
      hList = hMap->map[i];
      for (j=0, hNode=hList->first; j<hList->size; j++, hNode=hNode->next) {
	h = hNode->h;
	fprintf(file, "  %9.4f  %9.4f  %7.5f  %8.3f  %9.3e\n", h->pos[0]*RADIAN_TO_DEGREE, h->pos[1]*RADIAN_TO_DEGREE, h->z, h->w, h->M);
      }
    }
  }
  return;
}

void outAsciiFieldInfo(FILE *file, peak_param *pkPar)
{
  if (pkPar->field == 0) {
    fprintf(file, "# Field = %s, Omega = (%g, %g) [arcmin], theta_pix = %g [arcmin]\n", 
	    printField(pkPar->field), pkPar->Omega[0]*RADIAN_TO_ARCMIN, pkPar->Omega[1]*RADIAN_TO_ARCMIN, pkPar->theta_pix*RADIAN_TO_ARCMIN);
  }
  else if (pkPar->field == 1) {
    fprintf(file, "# Field = %s, Omega = (%g, %g) [arcmin], theta_pix = %g [arcmin]\n", 
	    printField(pkPar->field), pkPar->Omega[0]*RADIAN_TO_ARCMIN, pkPar->Omega[1]*RADIAN_TO_ARCMIN, pkPar->theta_pix*RADIAN_TO_ARCMIN);
    fprintf(file, "# nside = %ld, patch = %ld, rotAng = %g [deg], HP_resol = %d\n", pkPar->nside, pkPar->patch, pkPar->rotAng*RADIAN_TO_DEGREE, pkPar->HP_resol);
  }
  else {
    fprintf(file, "# Field = %s, nside = %ld, patch = %ld, HP_resol = %d\n", printField(pkPar->field), pkPar->nside, pkPar->patch, pkPar->HP_resol);
  }
  return;
}

void outAsciiHaloInfo(FILE *file, peak_param *pkPar)
{
  fprintf(file, "# z_halo_min = %g, z_halo_max = %g, N_z_halo = %d\n", pkPar->z_halo_min, pkPar->z_halo_max, pkPar->N_z_halo);
  fprintf(file, "# M_min = %8.2e, M_max = %8.2e [M_sol/h], dlogM = %g\n", pkPar->M_min, pkPar->M_max, pkPar->dlogM);
  return;
}

void outAsciiHaloCat(char name[], cosmo_hm *chPar, peak_param *pkPar, halo_map *hMap, error **err)
{
  if (pkPar->outHaloCat == 0) return;
  
  FILE *file = fopen_err(name, "w", err); forwardError(*err, __LINE__,);
  
  fprintf(file, "# Halo catalogue\n");
  fprintf(file, "#\n");
  outAsciiCosmoParam(file, chPar, pkPar);
  fprintf(file, "#\n");
  fprintf(file, "# Halo model parameters = %s\n", pkPar->hmParPath);
  fprintf(file, "# Mass function model = %s\n", smassfct_t(chPar->massfct));
  fprintf(file, "#\n");
  fprintf(file, "# Peak parameters = %s\n", pkPar->pkParPath);
  outAsciiFieldInfo(file, pkPar);
  outAsciiHaloInfo(file, pkPar);
  fprintf(file, "#\n");
  
  outAscii_halo_map(file, pkPar, hMap);
  
  fclose(file);
  if (pkPar->verbose < 3) printf("Outputed \"%s\"\n", name);
  return;
}

#ifdef __CAMELUS_USE_FITS__
void outFits_halo_t(FITS_t *fits, halo_t *h, double factor)
{
  float fBuff;
  fBuff = (float)(h->pos[0] * factor); writeTableColumn(fits, 0, 1, &fBuff);
  fBuff = (float)(h->pos[1] * factor); writeTableColumn(fits, 1, 1, &fBuff);
  fBuff = (float)h->z;                 writeTableColumn(fits, 2, 1, &fBuff);
  fBuff = (float)h->w;                 writeTableColumn(fits, 3, 1, &fBuff);
  fBuff = (float)h->M;                 writeTableColumn(fits, 4, 1, &fBuff);
  nextRow(fits);
  return;
}

void outFits_halo_map(FITS_t *fits, peak_param *pkPar, halo_map *hMap)
{
  if (pkPar->field == 0) {
    addColumn(fits, "THX", TFLOAT, "arcmin");
    addColumn(fits, "THY", TFLOAT, "arcmin");
  }
  else {
    addColumn(fits, "RA",  TFLOAT, "deg");
    addColumn(fits, "DEC", TFLOAT, "deg");
  }
  addColumn(fits, "Z", TFLOAT, "-       ");
  addColumn(fits, "W", TFLOAT, "Mpc/h");
  addColumn(fits, "M", TFLOAT, "M_sol/h");
  
  updateComment(fits, "TTYPE3", "Redshift");
  updateComment(fits, "TTYPE4", "Comving radial distance");
  updateComment(fits, "TTYPE5", "Halo mass");
  
  addLineSpread(fits);
  addKeyword(fits, TINT, "NBHALOS", &hMap->total, "[-] Number of halos");
  
  double factor = (pkPar->field == 0) ? RADIAN_TO_ARCMIN : RADIAN_TO_DEGREE;
  halo_list *hList;
  halo_node *hNode;
  int i, j;
  
  for (i=0; i<hMap->length; i++) {
    hList = hMap->map[i];
    for (j=0, hNode=hList->first; j<hList->size; j++, hNode=hNode->next) outFits_halo_t(fits, hNode->h, factor);
  }
  return;
}

void outFitsFieldInfo(FITS_t *fits, peak_param *pkPar)
{
  addLineSpread(fits);
  addKeyword(fits, TINT,    "FIELD",   &pkPar->field,    "[-] 0 = rectangle, 1 = projected HEALPix, 2 = HEALPix");
  
  double lfBuff;
  
  if (pkPar->field == 0) {
    lfBuff = pkPar->Omega[0] * RADIAN_TO_ARCMIN;
    addKeyword(fits, TDOUBLE, "OMEGAX",  &lfBuff,          "[arcmin] Field width");
    lfBuff = pkPar->Omega[1] * RADIAN_TO_ARCMIN;
    addKeyword(fits, TDOUBLE, "OMEGAY",  &lfBuff,          "[arcmin] Field height");
    lfBuff = pkPar->theta_pix * RADIAN_TO_ARCMIN;
    addKeyword(fits, TDOUBLE, "THPIX",   &lfBuff,          "[arcmin] Pixel size");
  }
  else if (pkPar->field == 1) {
    lfBuff = pkPar->Omega[0] * RADIAN_TO_ARCMIN;
    addKeyword(fits, TDOUBLE, "OMEGAX",  &lfBuff,          "[arcmin] Field width");
    lfBuff = pkPar->Omega[1] * RADIAN_TO_ARCMIN;
    addKeyword(fits, TDOUBLE, "OMEGAY",  &lfBuff,          "[arcmin] Field height");
    lfBuff = pkPar->theta_pix * RADIAN_TO_ARCMIN;
    addKeyword(fits, TDOUBLE, "THPIX",   &lfBuff,          "[arcmin] Pixel size");
    
    addKeyword(fits, TLONG,   "NSIDE",   &pkPar->nside,    "[-] N_side of the field");
    addKeyword(fits, TLONG,   "PATCH",   &pkPar->patch,    "[-] Ring order of the field");
    lfBuff = pkPar->rotAng * RADIAN_TO_DEGREE;
    addKeyword(fits, TDOUBLE, "ROTANG",  &lfBuff,          "[-] Rotation angle of the field after projection");
    addKeyword(fits, TINT,    "HPRESOL", &pkPar->HP_resol, "[-] Resolution for HEALPix regular galaxies or maps");
  }
  else {
    addKeyword(fits, TLONG,   "NSIDE",   &pkPar->nside,    "[-] N_side of the field");
    addKeyword(fits, TLONG,   "PATCH",   &pkPar->patch,    "[-] Ring order of the field");
    addKeyword(fits, TINT,    "HPRESOL", &pkPar->HP_resol, "[-] Resolution for HEALPix regular galaxies or maps");
  }
  return;
}

void outFitsHaloInfo(FITS_t *fits, peak_param *pkPar)
{
  char name[STRING_LENGTH_MAX];
  addLineSpread(fits);
  addKeyword(fits, TDOUBLE, "ZHALOMIN", &pkPar->z_halo_min, "[-] Minimum halo redshift");
  addKeyword(fits, TDOUBLE, "ZHALOMAX", &pkPar->z_halo_max, "[-] Maximum halo redshift");
  addKeyword(fits, TINT,    "NZHALO",   &pkPar->N_z_halo,   "[-] Number of halo redshift bins");
  sprintf(name, "%8.2e", pkPar->M_min);
  addKeyword(fits, TSTRING, "MMIN",     &name,              "[M_sol/h] Minimum halo mass");
  sprintf(name, "%8.2e", pkPar->M_max);
  addKeyword(fits, TSTRING, "MMAX",     &name,              "[M_sol/h] Maximum halo mass");
  addKeyword(fits, TDOUBLE, "DLOGM",    &pkPar->dlogM,      "[-] Halo mass binwidth");
  return;
}
#endif

void outFitsHaloCat(char name[], cosmo_hm *chPar, peak_param *pkPar, halo_map *hMap)
{
  if (pkPar->outHaloCat == 0) return;
  
#ifdef __CAMELUS_USE_FITS__
  FITS_t *fits = initializeTableWriter_FITS_t(name);
  outFits_halo_map(fits, pkPar, hMap);
  
  outFitsCosmoParam(fits, chPar, pkPar);
  
  addLineSpread(fits);
  addKeyword(fits, TSTRING, "HMPARAM", pkPar->hmParPath,           "[-] Path of the halo model parameter file");
  addKeyword(fits, TSTRING, "MASSFCT", smassfct_t(chPar->massfct), "[-] Mass function model");
  
  addLineSpread(fits);
  addKeyword(fits, TSTRING, "PKPARAM", pkPar->pkParPath,           "[-] Path of the peak parameter file");
  outFitsFieldInfo(fits, pkPar);
  outFitsHaloInfo(fits, pkPar);
  
  free_FITS_t(fits);
  if (pkPar->verbose < 3) printf("Outputed \"%s\"\n", name);
#endif
  return;
}

//----------------------------------------------------------------------
//-- Functions related to fast simulations

double dVol(cosmo_hm *chPar, double z, double ww, double dz, double area, error **err)
{
  //-- Comoving volume at z for a slice of dz and pkPar->area in [arcmin^2]
  //-- Since nicaea returns the mass function in comoving volume, here the comoving lightcone is used.
  
  double fK = f_K(chPar->cosmo, ww, err);                   forwardError(*err, __LINE__, -1.0);
  double dV;
  
  dV  = HUBBLE_DISTANCE * SQ(fK);
  dV *= area; //-- [rad^2]
  dV *= dz / sqrt(Esqr(chPar->cosmo, 1.0/(1.0+z), 0, err)); forwardError(*err, __LINE__, -1.0); //-- wOmegar = 0
  return dV;
}

void randomizePosition(peak_param *pkPar, double pos[2], error **err)
{
  double *Omega      = pkPar->Omega;
  gsl_rng *generator = pkPar->generator;
  
  //-- Set random position n in solid angle Omega
  if (pkPar->field == 0) {
    pos[0] = gsl_ran_flat(generator, 0.0, 1.0) * Omega[0];
    pos[1] = gsl_ran_flat(generator, 0.0, 1.0) * Omega[1];
  }
  else {
#ifdef __CAMELUS_USE_HEALPIX__
    patchSampling4(generator, pkPar->nside, pkPar->cap, pkPar->level, pkPar->length, pkPar->off, pkPar->j, pkPar->z0, pkPar->center[0], &pos[0], &pos[1]);
#endif
  }
  return;
}

void setMassSamplers(cosmo_hm *chPar, peak_param *pkPar, sampler_arr *hSampArr, int doVolume, error **err)
{
  if (pkPar->verbose < 2) {printf("Computing mass function...\r"); fflush(stdout);}
  
  //-- Structure for mass function
  cosmo_hm_params *cANDp = (cosmo_hm_params*)malloc_err(sizeof(cosmo_hm_params), err); forwardError(*err, __LINE__,);
  cANDp->model = chPar;
  cANDp->asymptotic = 0;
  
  int N_M   = pkPar->N_M;
  double dz = pkPar->dz_halo; 
  double r  = pow10(pkPar->dlogM); //-- Convert bin width to ratio
  sampler_t *hSamp;
  double *x, z, ww, dV;
  int i, j;
  
  //-- Loop over redshifts
  for (i=0, z=pkPar->z_halo_min+0.5*dz; i<hSampArr->length; i++, z+=dz) {
    hSamp     = hSampArr->array[i];
    hSamp->dx = pkPar->dlogM;
    x         = hSamp->x;
    x[N_M]    = pkPar->M_min;
    for (j=N_M; j>0; j--) x[j-1] = x[j] * r; //-- Fill masses in the reversed order to increase the sampling precision
    
    //-- Fill mass function and cdf
    cANDp->a = 1.0 / (1.0 + z);
    fillMassFct(cANDp, hSamp, err);                  forwardError(*err, __LINE__,);
    
    //-- Compute the volume
    if (doVolume == 1) { //-- This option is for kappa_1 calculation.
      ww = w(chPar->cosmo, cANDp->a, 0, err);        forwardError(*err, __LINE__,); //-- wOmegar = 0
      dV = dVol(chPar, z, ww, dz, pkPar->area, err); forwardError(*err, __LINE__,);
      hSamp->totPdf *= dV; //-- totPdf becomes the total number of halos in the volume dV
      hSamp->x_mean *= dV; //-- x_mean becomes the total mass in the volume dV
    }
  }
  
  free(cANDp);
  if (pkPar->verbose < 3) printf("Computed mass function    \n");
  return;
}

void rescaleWithVolume(cosmo_hm *chPar, peak_param *pkPar, sampler_arr *hSampArr, error **err)
{
  double dz = pkPar->dz_halo; 
  sampler_t *hSamp;
  double z, ww, dV;
  int i;
  
  //-- Rescale hSampArr with volume
  for (i=0, z=pkPar->z_halo_min+0.5*dz; i<hSampArr->length; i++, z+=dz) {
    hSamp = hSampArr->array[i];
    ww    = w(chPar->cosmo, 1.0/(1.0+z), 0, err);                 forwardError(*err, __LINE__,); //-- wOmegar = 0
    dV    = dVol(chPar, z, ww, pkPar->dz_halo, pkPar->area, err); forwardError(*err, __LINE__,);
    hSamp->totPdf *= dV; //-- totPdf becomes the total number of halos in the volume dV
    hSamp->x_mean *= dV; //-- x_mean becomes the total mass in the volume dV
  }
  return;
}

void sampleHalos(cosmo_hm *chPar, peak_param *pkPar, sampler_t *hSamp, halo_map *hMap, double z1, double z2, error **err)
{
  gsl_rng *generator = pkPar->generator;
  double z, p, M;
  
  //-- Compute the total number of halos and sample position
  int nbHalos = (int)round(hSamp->totPdf);
  double pos[2];
  
  //-- Sample masses and control
  int i, buff;
  for (i=0; i<nbHalos; i++) {
    p = gsl_ran_flat(generator, 0.0, 1.0);
    M = execute_sampler_t(hSamp, p);
    z = gsl_ran_flat(generator, z1, z2);
    randomizePosition(pkPar, pos, err);                        forwardError(*err, __LINE__,);
    append_halo_map(chPar, pkPar, hMap, pos, z, -1.0, M, err); forwardError(*err, __LINE__,);
  }
  return;
}

void makeFastSimul(cosmo_hm *chPar, peak_param *pkPar, sampler_arr *hSampArr, halo_map *hMap, error **err)
{
  //-- Reset
  reset_halo_map(hMap);
  
  double dz = pkPar->dz_halo;
  double z;
  int i;
  
  for (i=0, z=pkPar->z_halo_min; i<pkPar->N_z_halo; i++, z+=dz) {
    sampleHalos(chPar, pkPar, hSampArr->array[i], hMap, z, z+dz, err);
    forwardError(*err, __LINE__,);
  }
  
  if      (pkPar->verbose < 3)   printf("Generated %d halos         \n", hMap->total);
  else if (pkPar->verbose == 3) {printf("%6d halos, ", hMap->total); fflush(stdout);}
  return;
}

void readCatOrMakeSimulAndOutput(cosmo_hm *chPar, peak_param *pkPar, sampler_arr *hSampArr, halo_map *hMap, error **err)
{
  if (pkPar->doInHaloCat == 0) {
    makeFastSimul(chPar, pkPar, hSampArr, hMap, err);
    forwardError(*err, __LINE__,);
  }
  else {
    read_halo_map(pkPar->inHaloCatPath, chPar, pkPar, hMap, err);
    forwardError(*err, __LINE__,);
  }
  
  //-- Output
  char name[STRING_LENGTH_MAX];
  if (pkPar->doFITS == 0) {
    sprintf(name, "%shaloCat", pkPar->prefix);
    outAsciiHaloCat(name, chPar, pkPar, hMap, err); forwardError(*err, __LINE__,);
  }
  else {
    sprintf(name, "%shaloCat.fits", pkPar->prefix);
    outFitsHaloCat(name, chPar, pkPar, hMap);
  }
  return;
}

//----------------------------------------------------------------------
//-- Functions related to mass sheet

void setLambda(cosmo_hm *chPar, peak_param *pkPar, sampler_arr *hSampArr, double_arr *lambda)
{
  double rho_bar = CRITICAL_DENSITY * chPar->cosmo->Omega_m;
  int i;
  for (i=0; i<hSampArr->length; i++) lambda->array[i] = hSampArr->array[i]->x_mean / rho_bar;
  return;
}

void massSheet(cosmo_hm *chPar, peak_param *pkPar, sampler_t *hSamp, double sheet[2], error **err)
{
  //-- Compute the mass sheet values (kappa_0 & kappa_1) at a source redshift z_s
  
  if (pkPar->verbose < 2) {printf("Computing mass sheet...\r"); fflush(stdout);}
  
  int N_M  = pkPar->N_M;
  double r = pow(10, pkPar->dlogM); //-- Convert bin width to ratio
  double *x;
  int j;
  
  //-- Fill mass
  hSamp->dx = pkPar->dlogM;
  x         = hSamp->x;
  x[N_M]    = pkPar->M_min;
  for (j=N_M; j>0; j--) x[j-1] = x[j] * r; //-- Fill masses in the reversed order to increase the sampling precision
  
  //-- Structure for mass function
  cosmo_hm_params *cANDp = (cosmo_hm_params*)malloc_err(sizeof(cosmo_hm_params), err); forwardError(*err, __LINE__,);
  cANDp->model = chPar;
  cANDp->asymptotic = 0;
  
  double z_halo_min = pkPar->z_halo_min;
  double dz_halo    = pkPar->dz_halo;
  double z_s        = pkPar->z_s;
  double w_s        = w(chPar->cosmo, 1.0/(1.0+z_s), 0, err); forwardError(*err, __LINE__,); //-- wOmegar = 0, pkPar->w_s is not necessarily precomputed
  double kappa_0    = 0.0;
  double kappa_1    = 0.0;
  
  double fraction, value;
  double z_l, a_l, w_l;
  int i, i_max;
  
  i_max = (int)ceil((z_s - z_halo_min) / dz_halo);
  i_max = MIN(i_max, pkPar->N_z_halo);
  
  for (i=0; i<i_max; i++) {
    fraction = MIN(1.0, (z_s - z_halo_min) / dz_halo - i); //-- Deal with the fractional contribution in the last slice
    z_l      = z_halo_min + (i+fraction*0.5)*dz_halo;
    a_l      = 1.0 / (1.0 + z_l);
    
    //-- Fill mass function and cdf
    cANDp->a = a_l;
    fillMassFct(cANDp, hSamp, err);                                             forwardError(*err, __LINE__,);
    
    w_l      = w(chPar->cosmo, a_l, 0, err);                                    forwardError(*err, __LINE__,); //-- wOmegar = 0
    value    = f_K(chPar->cosmo, w_s - w_l, err) * f_K(chPar->cosmo, w_l, err); forwardError(*err, __LINE__,);
    value   *= fraction / (a_l * sqrt(Esqr(chPar->cosmo, a_l, 0, err)));        forwardError(*err, __LINE__,); //-- wOmegar = 0
    kappa_0 += value * CRITICAL_DENSITY * chPar->cosmo->Omega_m;
    kappa_1 += value * hSamp->x_mean; //-- x_mean = average mass density (comoving) contributed by all halos = lambda * rho_bar
  }
  free(cANDp);
  
  double factor = FOUR_PI_G_OVER_C2 * HUBBLE_DISTANCE * dz_halo / f_K(chPar->cosmo, w_s, err); forwardError(*err, __LINE__,);
  sheet[0] = kappa_0 * factor;
  sheet[1] = kappa_1 * factor;
  if (pkPar->verbose < 3) printf("Computed mass sheet    \n");
  return;
}

void setKappa1Interpolator(cosmo_hm *chPar, peak_param *pkPar, sampler_arr *hSampArr, interpolator_t *k1Inter, double *k0Arr, error **err)
{
  //-- hSampArr should already be set with doVolume = 0
  
  int N_z_gal       = k1Inter->length -1;
  double z_halo_min = pkPar->z_halo_min;
  double z_gal_min  = chPar->redshift->par_nz[0];
  double dz_halo    = pkPar->dz_halo;
  double dz_gal     = pkPar->dz_gal;
  double factor     = FOUR_PI_G_OVER_C2 * HUBBLE_DISTANCE * dz_halo;
  cosmo *cmPar      = chPar->cosmo;
  
  k1Inter->dx = dz_gal;
  
  double kappa_0, kappa_1, fraction, kernel;
  double z_s, a_s, w_s, fK_inv_s;
  double z_l, a_l, w_l;
  int i, j, i_max;
  
  for (j=0, z_s=z_gal_min; j<=N_z_gal; j++, z_s+=dz_gal) {
    kappa_0 = 0.0;
    kappa_1 = 0.0;
    
    k1Inter->x[j] = z_s;
    a_s  = 1.0 / (1.0 + z_s);
    w_s  = w(cmPar, a_s, 0, err);          forwardError(*err, __LINE__,); //-- wOmegar = 0
    fK_inv_s = 1.0 / f_K(cmPar, w_s, err); forwardError(*err, __LINE__,);
    
    i_max = (int)ceil((z_s - z_halo_min) / dz_halo);
    i_max = MIN(i_max, pkPar->N_z_halo);
    
    for (i=0; i<i_max; i++) {
      fraction = MIN(1.0, (z_s - z_halo_min) / dz_halo - i); //-- Deal with the fractional contribution in the last slice
      z_l      = z_halo_min + (i+fraction*0.5)*dz_halo;
      a_l      = 1.0 / (1.0 + z_l);
      w_l      = w(cmPar, a_l, 0, err);                                           forwardError(*err, __LINE__,); //-- wOmegar = 0
      kernel   = f_K(cmPar, w_s - w_l, err) * f_K(cmPar, w_l, err) * (1.0 + z_l); forwardError(*err, __LINE__,);
      kernel  *= fraction / sqrt(Esqr(cmPar, a_l, 0, err));                       forwardError(*err, __LINE__,); //-- wOmegar = 0
      kappa_0 += kernel * CRITICAL_DENSITY * cmPar->Omega_m; //-- kernel * rho_bar
      kappa_1 += kernel * hSampArr->array[i]->x_mean;        //-- kernel * rho_bar * lambda
    }
    
    if (k0Arr != NULL) k0Arr[j] = kappa_0 * factor * fK_inv_s;
    k1Inter->value[j] = kappa_1 * factor * fK_inv_s;
  }
  return;
}

void setMassSampAndK1Inter(cosmo_hm *chPar, peak_param *pkPar, sampler_arr *hSampArr, double_arr *lambda, interpolator_t *k1Inter, error **err)
{
  if (pkPar->doSubtraction == 2) {
    setMassSamplers(chPar, pkPar, hSampArr, 0, err);                   forwardError(*err, __LINE__,); //-- doVolume = 0
    if (lambda != NULL) setLambda(chPar, pkPar, hSampArr, lambda);
    setKappa1Interpolator(chPar, pkPar, hSampArr, k1Inter, NULL, err); forwardError(*err, __LINE__,); //-- k0Arr = NULL
    rescaleWithVolume(chPar, pkPar, hSampArr, err);                    forwardError(*err, __LINE__,);
  }
  else if (pkPar->doInHaloCat == 0) {
    setMassSamplers(chPar, pkPar, hSampArr, 1, err);                   forwardError(*err, __LINE__,); //-- doVolume = 1
  }
  return;
}

void outAsciiLambda(char name[], cosmo_hm *chPar, peak_param *pkPar, double_arr *lambda, error **err)
{
  FILE *file = fopen_err(name, "w", err); forwardError(*err, __LINE__,);
  
  fprintf(file, "# lambda - proportion of total mass contributed by halos\n");
  fprintf(file, "#\n");
  outAsciiCosmoParam(file, chPar, pkPar);
  fprintf(file, "#\n");
  fprintf(file, "# Peak parameters = %s\n", pkPar->pkParPath);
  fprintf(file, "# z_halo_min = %f, z_halo_max = %f, N_z_halo = %d\n", pkPar->z_halo_min, pkPar->z_halo_max, pkPar->N_z_halo);
  fprintf(file, "# M_min = %8.2e, M_max = %8.2e [M_sol/h]\n", pkPar->M_min, pkPar->M_max);
  fprintf(file, "#\n");
  fprintf(file, "#      z        lambda\n");
  fprintf(file, "#     [-]          [-]\n");
  
  double dz = pkPar->dz_halo;
  double z;
  int i;
  
  for (i=0, z=pkPar->z_halo_min+0.5*dz; i<lambda->length; i++, z+=dz) {
    fprintf(file, "  %7.5f  %11.5e\n", z, lambda->array[i]);
  }
  
  fclose(file);
  if (pkPar->verbose < 3) printf("Outputed \"%s\"\n", name);
  return;
}

void outAsciiMassSheet(char name[], cosmo_hm *chPar, peak_param *pkPar, interpolator_t *k1Inter, double *k0Arr, error **err)
{
  FILE *file = fopen_err(name, "w", err); forwardError(*err, __LINE__,);
  
  fprintf(file, "# Mass sheets\n");
  fprintf(file, "#\n");
  outAsciiCosmoParam(file, chPar, pkPar);
  fprintf(file, "#\n");
  fprintf(file, "# Halo model parameters = %s\n", pkPar->hmParPath);
  fprintf(file, "# z_s_min = %f, z_s_max = %f\n", chPar->redshift->par_nz[0], chPar->redshift->par_nz[1]);
  fprintf(file, "#\n");
  fprintf(file, "# Peak parameters = %s\n", pkPar->pkParPath);
  fprintf(file, "# z_halo_min = %f, z_halo_max = %f, N_z_halo = %d\n", pkPar->z_halo_min, pkPar->z_halo_max, pkPar->N_z_halo);
  fprintf(file, "# M_min = %8.2e, M_max = %8.2e [M_sol/h]\n", pkPar->M_min, pkPar->M_max);
  fprintf(file, "#\n");
  fprintf(file, "#      z       kappa_0      kappa_1\n");
  fprintf(file, "#     [-]          [-]          [-]\n");
  
  int i;
  for (i=0; i<k1Inter->length; i++) {
    fprintf(file, "  %7.5f  %11.5e  %11.5e\n", k1Inter->x[i], k0Arr[i], k1Inter->value[i]);
  }
  
  fclose(file);
  if (pkPar->verbose < 3) printf("Outputed \"%s\"\n", name);
  return;
}

//----------------------------------------------------------------------
//-- Main functions

void doMassFct(cosmo_hm *chPar, peak_param *pkPar, double z, error **err)
{
  //-- Default values
  int N = 10;
  double z_min = 0.05;
  double dz = 0.1;
  
  char name[STRING_LENGTH_MAX];
  char buffer[STRING_LENGTH_MAX];
  double y;
  int i;
  
  if (z < 0.0) {
    printf("Want to generate all files at once. Lazy, hein? I understand.\n");
    printf("Enter the lowest redshift (default: 0.05): ");
    fgets(buffer, STRING_LENGTH_MAX, stdin);
    if (strcmp(buffer, "\n")) z_min = atof(buffer);
    
    printf("Enter dz (default: 0.1): ");
    fgets(buffer, STRING_LENGTH_MAX, stdin);
    if (strcmp(buffer, "\n")) dz = atof(buffer);
    
    printf("Enter the number of outputs (default: 10): ");
    fgets(buffer, STRING_LENGTH_MAX, stdin);
    if (strcmp(buffer, "\n")) N = atoi(buffer);
    
    printf("Enter the prefix of outputs (default: camelus_massFct_z): "); 
    fgets(buffer, STRING_LENGTH_MAX, stdin);
    if (strcmp(buffer, "\n")) strtok(buffer, "\n");
    else sprintf(buffer, "camelus_massFct_z");
    
    printf("\nComputing mass function...\r"); fflush(stdout);
    for (i=0, y=z_min; i<N; i++, y+=dz) {
      if (pkPar->doFITS == 0) {
	sprintf(name, "%s%g", buffer, y);
	outAsciiMassFct(name, chPar, pkPar, y, err); forwardError(*err, __LINE__,);
      }
      else {
	sprintf(name, "%s%g.fits", buffer, y);
	outFitsMassFct(name, chPar, pkPar, y, err); forwardError(*err, __LINE__,);
      }
    }
  }
  
  else {
    if (pkPar->doFITS == 0) {
      sprintf(name, "%smassFct_%s_z%.3f", pkPar->prefix, smassfct_t(chPar->massfct), z);
      outAsciiMassFct(name, chPar, pkPar, z, err);  forwardError(*err, __LINE__,);
    }
    else {
      sprintf(name, "%smassFct_%s_z%.3f.fits", pkPar->prefix, smassfct_t(chPar->massfct), z);
      outFitsMassFct(name, chPar, pkPar, z, err);  forwardError(*err, __LINE__,);
    }
  }
  
  printf("------------------------------------------------------------------------\n");
  return;
}

void doMassSheet(cosmo_hm *chPar, peak_param *pkPar, double z_s, error **err)
{
  sampler_t *hSamp;
  sampler_arr *hSampArr;
  double_arr *lambda;
  interpolator_t *k1Inter;
  char name[STRING_LENGTH_MAX];
  double k0Arr[pkPar->N_z_gal+1];
  double sheet[2];
  
  //-- Compute mass sheet at a redshift
  if (z_s > 0.0) {
    pkPar->z_s = z_s;
    
    printf("Source at z_s = %f\n", pkPar->z_s);
    printf("z_halo_min = %f, z_halo_max = %f, %d bins of dz = %f\n", pkPar->z_halo_min, pkPar->z_halo_max, pkPar->N_z_halo, pkPar->dz_halo);
    printf("M_min = %8.2e, M_max = %8.2e [M_sol/h]\n", pkPar->M_min, pkPar->M_max);
    printf("\n");
    
    hSamp = initialize_sampler_t(pkPar->N_M+1);
    massSheet(chPar, pkPar, hSamp, sheet, err); forwardError(*err, __LINE__,);
    
    printf("For the given configuration:\n");
    printf("kappa_0 = %f\n", sheet[0]);
    printf("kappa_1 = %f\n", sheet[1]);
    //printf("sigma_noise = %f\n", pkPar->sigma_noise[0]);
    free_sampler_t(hSamp);
  }
  
  //-- Compute mass sheet for a series of redshifts
  else {
    hSampArr = initialize_sampler_arr(pkPar->N_z_halo, pkPar->N_M+1);
    lambda   = initialize_double_arr(pkPar->N_z_halo);
    k1Inter  = initialize_interpolator_t(pkPar->N_z_gal+1);
    
    setMassSamplers(chPar, pkPar, hSampArr, 0, err);                    forwardError(*err, __LINE__,); //-- doVolume = 0
    setLambda(chPar, pkPar, hSampArr, lambda);
    setKappa1Interpolator(chPar, pkPar, hSampArr, k1Inter, k0Arr, err); forwardError(*err, __LINE__,);
    
    sprintf(name, "%slambda", pkPar->prefix);
    outAsciiLambda(name, chPar, pkPar, lambda, err);                      forwardError(*err, __LINE__,);
    sprintf(name, "%smassSheet", pkPar->prefix);
    outAsciiMassSheet(name, chPar, pkPar, k1Inter, k0Arr, err);           forwardError(*err, __LINE__,);
    
    free_sampler_arr(hSampArr);
    free_double_arr(lambda);
    free_interpolator_t(k1Inter);
  }
  
  printf("------------------------------------------------------------------------\n");
  return;
}

void doFastSimulation(cosmo_hm *chPar, peak_param *pkPar, error **err)
{
  sampler_arr *hSampArr = initialize_sampler_arr(pkPar->N_z_halo, pkPar->N_M+1);
  setMassSamplers(chPar, pkPar, hSampArr, 1, err);                                                      forwardError(*err, __LINE__,); //-- doVolume = 1
  halo_map *hMap        = initialize_halo_map(pkPar->resol[0], pkPar->resol[1], pkPar->theta_pix, err); forwardError(*err, __LINE__,);
  
  readCatOrMakeSimulAndOutput(chPar, pkPar, hSampArr, hMap, err); forwardError(*err, __LINE__,);
  
  free_sampler_arr(hSampArr);
  free_halo_map(hMap);
  printf("------------------------------------------------------------------------\n");
  return;
}

//----------------------------------------------------------------------

