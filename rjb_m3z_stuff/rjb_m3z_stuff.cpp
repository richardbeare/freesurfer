/**
 * @file  rjb_m3z_stuff.cpp
 * @brief A program to explore warp files
 *
 */

/* Designed to convert m3z files derived from fnirt/itk warps to 
   a form appropriate for use my mri_ca_register.

   Uses the same technique as mri_concatenate_warp to appropriately
   set source and target geometeries.

   Generates an affine matrix, from the talairach.lta file to attach.
   This appears important for mri_ca_register

   Resamples the gca nodes to a 256x256x256 grid. Seems the nodes need
   to match the target size.

*/

/*
 * Original Author: Richard Beare
 *
 * mri_convert_warp used as template
 */

#include <string>
#include <iostream>
#include <fstream>
#include <iomanip>
// C compat
#ifdef __cplusplus
extern "C"
{
#endif

#include "error.h"
#include "gcamorph.h"
#include "macros.h"
#include "matrix.h"
#include "mri.h"
#include "mri_circulars.h"
#include "version.h"
#include "diag.h"
#include "cma.h"

#ifdef __cplusplus
}
#endif

using namespace std;
namespace filetypes {
enum FileType { UNKNOWN, M3Z, FSL, ITK };
}

struct Parameters
{
  string in_warp;
  string out_warp;
  string affineMat;
  string sourceImage;
  string atlasImage;
  string maskImage;
  filetypes::FileType in_type;
  filetypes::FileType out_type;
  bool infoOnly;
};

static struct Parameters P =
  { "", "", "", "", "", "", filetypes::UNKNOWN, filetypes::UNKNOWN, false};

static void printUsage(void);
static bool parseCommandLine(int argc, char *argv[], Parameters & P);

static char vcid[] =
    "$Id: rjb_m3z_stuff.cpp,v 1.1 2018/11/01 $";
const char *Progname = NULL;

GCAM* readM3Z(const string& warp_file)
// Read an m3z file. Just calls down to GCAMread
{
  GCAM* gcam = GCAMread(warp_file.c_str());
  if (gcam == NULL)
  {
    cerr << "ERROR readM3Z: cannot read " << warp_file << endl;
    exit(1);
  }

  return gcam;
}

/* from mri_info */
int PrettyMatrixPrint(const MATRIX *mat)
{

  if (mat == NULL)
  {
    ErrorReturn(ERROR_BADPARM,(ERROR_BADPARM, "mat = NULL!")) ;
  }

  if (mat->type != MATRIX_REAL)
  {
    ErrorReturn(ERROR_BADPARM,(ERROR_BADPARM, "mat is not MATRIX_REAL type")) ;
  }

  if (mat->rows != 4 || mat->cols != 4)
  {
    ErrorReturn(ERROR_BADPARM,(ERROR_BADPARM, "mat is not of 4 x 4")) ;
  }
#if 0
  for (row=1; row < 5; ++row)
  {
    printf("              %8.4f %8.4f %8.4f %10.4f\n",
           mat->rptr[row][1], mat->rptr[row][2],
           mat->rptr[row][3], mat->rptr[row][4]);
  }
#endif
  MatrixPrint(stdout, mat);
  return (NO_ERROR);
}

void GCAMInitStuff(GCA_MORPH *gcam, GCA *gca, string MRIf="")
{
    // added by xhan
    int x, y, z, n, label, max_n, max_label;
    float max_p;
    GC1D *gc;
    GCA_MORPH_NODE  *gcamn ;
    GCA_PRIOR *gcap;

    //gcam->ninputs = mri_inputs->nframes ;
    //getVolGeom(mri_inputs, &gcam->image);
    GCAsetVolGeom(gca, &gcam->atlas);
    gcam->gca = gca ;
    gcam->spacing = gca->prior_spacing;

    // use gca information
    for (x = 0 ; x < gcam->width ; x++)
    {
      for (y = 0 ; y < gcam->height ; y++)
      {
        for (z = 0 ; z < gcam->depth ; z++)
        {
          gcamn = &gcam->nodes[x][y][z] ;
          gcap = &gca->priors[x][y][z] ;
          max_p = 0 ;
          max_n = -1 ;
          max_label = 0 ;

          // find the label which has the max p
          for (n = 0 ; n < gcap->nlabels ; n++)
          {
            label = gcap->labels[n] ;   // get prior label
            if (label == Gdiag_no)
            {
              DiagBreak() ;
            }
            if (label >= MAX_CMA_LABEL)
            {
              printf("invalid label %d at (%d, %d, %d) in prior volume\n",
                     label, x, y, z);
            }
            if (gcap->priors[n] >= max_p) // update the max_p and max_label
            {
              max_n = n ;
              max_p = gcap->priors[n] ;
              max_label = gcap->labels[n] ;
            }
          }

          gcamn->label = max_label ;
          gcamn->n = max_n ;
          gcamn->prior = max_p ;
          gc = GCAfindPriorGC(gca, x, y, z, max_label) ;
          // gc can be NULL
          gcamn->gc = gc ;
          gcamn->log_p = 0 ;

        }
    }

    GCAMcomputeOriginalProperties(gcam) ;
    if (MRIf != "")
    {
      MRI *mri_inputs = MRIread(MRIf.c_str());
      GCAMcomputeLabels(mri_inputs, gcam) ;
      MRIfree(&mri_inputs);
    }
    else
    {
     GCAMcomputeMaxPriorLabels(gcam) ;
    }
  }
}
void dumpMRI(MRI *M)
{

  cout << "Value" << MRIgetVoxVal(M, 0, 0, 0, 0) << endl;
}

void VGInfo(const VOL_GEOM *vg)
{
  cout << "VG Valid: " << vg->valid << endl;
  cout << "width " << vg->width << ", height " << vg->height << ", depth " << vg->depth << endl;
  cout << "xsize " << vg->xsize << ", ysize " << vg->ysize << ", zsize " << vg->zsize << endl;
  cout << "fname " << vg->fname << endl;

  cout << std::fixed << std::setw( 11 ) << std::setprecision( 6 ) << 
    "xras " << vg->x_r << " " << vg->x_a << " " << vg->x_s << endl;
  cout << std::fixed << std::setw( 11 ) << std::setprecision( 6 ) <<
    "yras " << vg->y_r << " " << vg->y_a << " " << vg->y_s << endl;
  cout << std::fixed << std::setw( 11 ) << std::setprecision( 6 ) << 
    "zras " << vg->z_r << " " << vg->z_a << " " << vg->z_s << endl;
  cout << std::fixed << std::setw( 11 ) << std::setprecision( 6 ) << 
    "cras " << vg->c_r << " " << vg->c_a << " " << vg->c_s << endl;

}

MRI *conformedButNotCoronal(MRI *sampleforslicedirection)
{
  MRI *templ;
  templ = MRIallocHeader(256, 256, 256, MRI_UCHAR, 1);

  templ->imnr0 = 1;
  templ->imnr1 = 256;
  templ->thick = 1.0;
  templ->ps = 1.0;
  templ->xsize = templ->ysize = templ->zsize = 1.0;
  templ->xstart = templ->ystart = templ->zstart = -128.0;
  templ->xend = templ->yend = templ->zend = 128.0;

  int sliceDirection = getSliceDirection(sampleforslicedirection);
  setDirectionCosine(templ, sliceDirection);
  // templ->c_r = sampleforslicedirection->c_r;
  // templ->c_a = sampleforslicedirection->c_a;
  // templ->c_s = sampleforslicedirection->c_s;
  // these are the template settings we are after
  templ->c_r = 0;
  templ->c_a = 0;
  templ->c_s = 0;

  templ->ras_good_flag = 1;
  templ->tr = sampleforslicedirection->tr;
  templ->te = sampleforslicedirection->te;
  templ->flip_angle = sampleforslicedirection->flip_angle;
  templ->ti = sampleforslicedirection->ti;
  return(templ);
}
MRI *conformedCoronal()
{
  MRI *templ;
  templ = MRIallocHeader(256, 256, 256, MRI_UCHAR, 1);

  templ->imnr0 = 1;
  templ->imnr1 = 256;
  templ->thick = 1.0;
  templ->ps = 1.0;
  templ->xsize = templ->ysize = templ->zsize = 1.0;
  templ->xstart = templ->ystart = templ->zstart = -128.0;
  templ->xend = templ->yend = templ->zend = 128.0;

  // int sliceDirection = getSliceDirection(sampleforslicedirection);
  // setDirectionCosine(templ, sliceDirection);
  setDirectionCosine(templ, MRI_CORONAL);
  // templ->c_r = sampleforslicedirection->c_r;
  // templ->c_a = sampleforslicedirection->c_a;
  // templ->c_s = sampleforslicedirection->c_s;
  // these are the template settings we are after
  templ->c_r = 0;
  templ->c_a = 0;
  templ->c_s = 0;

  templ->ras_good_flag = 1;
  return(templ);
}

void dumpGCAM(GCAM *gcam) {

  if (gcam->m_affine == NULL) {
    cout << "WARN - affine matrix null - not dumping gcam" << endl;
    return;
  }

  const VOL_GEOM *atlas = &(gcam->atlas);
  VECTOR *v1, *v2, *v3;
  v1 = VectorAlloc(4, MATRIX_REAL);
  *MATRIX_RELT(v1, 4, 1) = 1.0;
  v2 = MatrixCopy(v1, NULL);
  v3 = MatrixCopy(v1, NULL);

  MATRIX *atlasvox2ras = MatrixInverse(gcam->m_affine, NULL);
  
  for(int s=0; s < gcam->depth; s++)
  {
    for(int c=0; c < gcam->width; c++)
    {
      for(int r=0; r < gcam->height; r++)
      {
	GCA_MORPH_NODE* node = &gcam->nodes[c][r][s];
	float px, py, pz;
	float x = c * gcam->spacing;
	float y = r * gcam->spacing;
	float z = s * gcam->spacing;
	V3_X(v1) = x;
	V3_Y(v1) = y;
	V3_Z(v1) = z;
	MatrixMultiply(atlasvox2ras, v1, v2);
	px = V3_X(v2);
	py = V3_Y(v2);
	pz = V3_Z(v2);

	//	if (node->origx != node->x) {
	printf("Node: (%d, %d, %d), n(%f, %f, %f), ox(%f, %f, %f), (%f, %f, %f)\n",
	       c, r, s,
	       px, py, pz,
	       node->origx, node->origy, node->origz,
	       node->x, node->y, node->z);
	//	}
      }
    }
  }

}

void Info(const GCAM* gcam)
{
  cout << "Width, height, depth " << gcam->width << " " << gcam->height << " " << gcam->depth << endl;
  cout << "GCA pointer: " << gcam->gca << endl;
  cout << "Node pointer: " << gcam->nodes << endl;
  cout << "Atlas VG" << endl;
  VGInfo(&(gcam->atlas));
  cout << endl;
  cout << "Image VG" << endl;
  VGInfo(&(gcam->image));
  PrettyMatrixPrint((gcam->m_affine));

  cout << "mri_xind: " << gcam->mri_xind << endl; 
  
  
  cout << "Type : " << gcam->type << endl;
  cout << "status : " << gcam->status << endl;

}

MATRIX *MRIgetConformMatrixPad(const VOL_GEOM *vg, const MRI *target)
{
  // matrix describing the padding to produce a 256^3 image, but without
  // setting to coronal direction cosines
  MRI *mri;
  MATRIX *m_resample;

  mri = MRIallocFromVolGeom(const_cast<VOL_GEOM *>(vg), MRI_UCHAR, 1, 1);
  m_resample = MRIgetResampleMatrix(mri, const_cast<MRI *>(target));
  
  return (m_resample);

}

void maskGCAM(GCAM * gcam, MRI *mask)
{
  const VOL_GEOM *atlas = &(gcam->atlas);
  VECTOR *v1, *v2;
  v1 = VectorAlloc(4, MATRIX_REAL);
  *MATRIX_RELT(v1, 4, 1) = 1.0;
  v2 = MatrixCopy(v1, NULL);

  MATRIX *atlasvox2ras = MatrixInverse(gcam->m_affine, NULL);
  for(int s=0; s < gcam->depth; s++)
  {
    for(int c=0; c < gcam->width; c++)
    {
      for(int r=0; r < gcam->height; r++)
      {
	float x = c * gcam->spacing;
	float y = r * gcam->spacing;
	float z = s * gcam->spacing;
 
        float mval = 1;
	mval = MRIgetVoxVal(mask, x, y, z, 0);

	if (mval == 0) {
	  GCA_MORPH_NODE* node = &gcam->nodes[c][r][s];
	  float px, py, pz;

	  V3_X(v1) = x;
	  V3_Y(v1) = y;
	  V3_Z(v1) = z;
	  MatrixMultiply(atlasvox2ras, v1, v2);
	  px = V3_X(v2);
	  py = V3_Y(v2);
	  pz = V3_Z(v2);

	  node->invalid = GCAM_AREA_INVALID;
          node->origx = px;
          node->origy = py;
          node->origz = pz;
          node->x = px;
          node->y = py;
          node->z = pz;
	}
      }
    }
  }

}
  

void resampleGCAMto256(const GCAM *fnirt_gcam, GCAM *result_gcam)
{
  // the fnirt gcam isn't 256^3, make it so
  // Also need to swap the dimensions in the gcamorph nodes around, if
  // the image axes have changes, as per conforming process.
  //
  // lps2ras conversion is aleady done in mri_warp_convert
  int r, c, s;
  const VOL_GEOM *fnirt_atlas = &(fnirt_gcam->atlas);
  VECTOR *v1, *v2, *v3;

  /* an MRI header for the fnirt atlas */
  // MRI *fnirt_atlas_header = MRIallocFromVolGeom(const_cast<VOL_GEOM *>(fnirt_atlas), MRI_UCHAR, 1, 1);
  /* an MRI header for the fnirt atlas in a 256^3 volume */
  // MRI *fnirt_atlas_new = conformedButNotCoronal(fnirt_atlas_header);
  MRI *fs_tmplt = conformedCoronal();
  /* vox2vox transform between the original fnirt atlas and the new one */
  //MATRIX *conformMat = MRIgetConformMatrixPad(fnirt_atlas, fnirt_atlas_new);
  //MATRIX *conformMat = MRIgetConformMatrixPad(fnirt_atlas, fs_tmplt);
  VOL_GEOM sortofconformed;
  getVolGeom(fs_tmplt, &sortofconformed);
  MATRIX *conformCoronalMat = MRIgetConformMatrixPad(fnirt_atlas, fs_tmplt);

  // make an axis swap matrix for entries in the gcamorph structure.
  // This ought to be the conform matrix without the signs or offsets.
  MATRIX *gcamswap = MatrixCopy(conformCoronalMat, NULL);

  for (unsigned row=1; row < 4; ++row)
    {
      gcamswap->rptr[row][1] = fabs(gcamswap->rptr[row][1]);
      gcamswap->rptr[row][2] = fabs(gcamswap->rptr[row][2]);
      gcamswap->rptr[row][3] = fabs(gcamswap->rptr[row][3]);
      gcamswap->rptr[row][4] = 0;
    }
  // not yet sure if I need the inverse of this.
  
  // need to adust the headers in the new atlas
  PrettyMatrixPrint(MRImatrixOfDirectionCosines(fs_tmplt, NULL));
  cout << "=======================" << endl;
  PrettyMatrixPrint(conformCoronalMat);
  cout << "=======================" << endl;
  PrettyMatrixPrint(gcamswap);

  
  cout << "Starting resample" << endl;


  
  int atlaswidth = fnirt_atlas->width;
  int atlasheight = fnirt_atlas->height;
  int atlasdepth = fnirt_atlas->depth;

  v1 = VectorAlloc(4, MATRIX_REAL);
  *MATRIX_RELT(v1, 4, 1) = 1.0;
  v2 = MatrixCopy(v1, NULL);
  v3 = MatrixCopy(v1, NULL);
  
  V3_X(v1) = 11;
  V3_Y(v1) = 23;
  V3_Z(v1) = 37;

  MatrixMultiply(conformCoronalMat, v1, v2);

  cout << "(" << V3_X(v2) << ", " << V3_Y(v2) << ", " << V3_Z(v2) << ")" << endl;

  cout << "Starting loop" << endl;

  for (c = 0; c < result_gcam->width; c++) {
    for (r = 0; r < result_gcam->height; r++) {
      for (s = 0; s < result_gcam->depth; s++) {
	// New coordinates of the node in the 256^3 volume
	V3_X(v1) = c;
        V3_Y(v1) = r;
        V3_Z(v1) = s;
	// where does it come from in the original node map
	MatrixMultiply(conformCoronalMat, v1, v2);
	int origc = V3_X(v2);
	int origr = V3_Y(v2);
	int origs = V3_Z(v2);

	GCA_MORPH_NODE *newnode = &(result_gcam->nodes[c][r][s]);
	if (origc < 0 || origc > (atlaswidth - 1) ||
	    origr < 0 || origr > (atlasheight - 1) ||
	    origs < 0 || origs > (atlasdepth - 1)) {
	  // an invalid node
	  newnode->origx = 0;
	  newnode->origy = 0;
	  newnode->origz = 0;

	  newnode->x = 0;
	  newnode->y = 0;
	  newnode->z = 0;
          newnode->invalid = GCAM_POSITION_INVALID;
	  
	} else {
	  // copy from original
	  GCA_MORPH_NODE *oldnode = &(fnirt_gcam->nodes[origc][origr][origs]);
	  newnode->xn = c;
	  newnode->yn = r;
	  newnode->zn = s;
	  // fprintf(stdout, "(%d, %d, %d) (%d, %d, %d)\n", c,r,s, origc, origr, origs);
	  // fprintf(stdout, "(%f, %f, %f)(%f, %f, %f)\n",
	  //   	  oldnode->x, oldnode->y, oldnode->z,
	  // 	  oldnode->origx,oldnode->origy,oldnode->origz);

#if 0
	  V3_X(v2) = oldnode->x;
	  V3_Y(v2) = oldnode->y;
	  V3_Z(v2) = oldnode->z;
	  MatrixMultiply(gcamswap, v2, v3);

	  newnode->x = V3_X(v3);
	  newnode->y = V3_Y(v3);
	  newnode->z = V3_Z(v3);
	  
	  V3_X(v2) = oldnode->origx;
	  V3_Y(v2) = oldnode->origy;
	  V3_Z(v2) = oldnode->origz;
	  MatrixMultiply(gcamswap, v2, v3);
	  
	  newnode->origx = V3_X(v3);
	  newnode->origy = V3_Y(v3);
	  newnode->origz = V3_Z(v3);
#else
	  newnode->x = oldnode->x;
	  newnode->y = oldnode->y;
	  newnode->z = oldnode->z;

	  newnode->origx = oldnode->origx;
	  newnode->origy = oldnode->origy;
	  newnode->origz = oldnode->origz;

#endif
	  // fprintf(stdout, "nn(%f, %f, %f)(%f, %f, %f)\n",
	  //   	  newnode->x, newnode->y, newnode->z,
	  // 	  newnode->origx, newnode->origy, newnode->origz);
	}
      }
    }
  }

  cout << "finished loop" << endl;

  strcpy(sortofconformed.fname, "dummy");

  copyVolGeom(const_cast<const VOL_GEOM *>(&sortofconformed), &result_gcam->atlas);
  copyVolGeom(&fnirt_gcam->image, &result_gcam->image);
  
  result_gcam->m_affine = fnirt_gcam->m_affine;
  result_gcam->det = fnirt_gcam->det;
  result_gcam->type = fnirt_gcam->type;
  result_gcam->spacing = fnirt_gcam->spacing;
}

#if 0
void resampleGCAM(const GCAM *fnirt_gcam, const GCAM *fnirtvolchange_gcam, GCAM *result_gcam)
{
  // fnirt_gcam has the incorrect FOV, but produces correct warp results with mri_convert
  // i.e the result aligns nicely with the atlas. However the FOV is wrong and this
  // causes problem. fnirtvolchange_gcam has been informed of the correct src and
  // volume images, but updating the morph nodes doesn't work properly because of the differing
  // dimensions.
  // result_gcam has been initialised to the correct size, but is empty.
  // For each voxel in the atlas, we need to create the appropriate node in result_gcam.
  // This means figuring out the world coordinates, then finding the equivalent node in fnirt_gcam.
  // This is specifically for between fsl and fs templates, so we don't need to worry about
  // interpolating values anywhere - there is voxel equivalence. Only the FOVs differ.
  int r, c, s;
  const VOL_GEOM *fs_atlas = &(fnirtvolchange_gcam->atlas);
  const VOL_GEOM *fnirt_atlas = &(fnirt_gcam->atlas);

  MATRIX *conformMat = MRIgetConformMatrixPad(fnirt_atlas);
  MATRIX *fs_vox2ras = vg_getVoxelToRasXform(fs_atlas);
  MATRIX *fnirt_ras2vox = vg_getRasToVoxelXform(fnirt_atlas);
  VECTOR *v1, *v2, *v3;

  PrettyMatrixPrint(fs_vox2ras);
  cout << "---------------------" << endl;
  PrettyMatrixPrint(vg_getVoxelToRasXform(fnirt_atlas));
  cout << "---------------------" << endl;
  PrettyMatrixPrint(conformMat);
  
  v1 = VectorAlloc(4, MATRIX_REAL);
  *MATRIX_RELT(v1, 4, 1) = 1.0;
  v2 = MatrixCopy(v1, NULL);
  v3 = MatrixCopy(v1, NULL);

  V3_X(v1) = 11;
  V3_Y(v1) = 23;
  V3_Z(v1) = 37;

  MatrixMultiply(fs_vox2ras, v1, v2);
  // v2 contains world coordinates
  MatrixMultiply(fnirt_ras2vox, v2, v3);

  cout << "(" << V3_X(v3) << ", " << V3_Y(v3) << ", " << V3_Z(v3) << ")" << endl;
  

  for (c = 0; c < fs_atlas->width; c++) {
    for (r = 0; r < fs_atlas->height; r++) {
      for (s = 0; s < fs_atlas->depth; s++) {
	// world coordinates of this voxel
	V3_X(v1) = c;
        V3_Y(v1) = r;
        V3_Z(v1) = s;

	MatrixMultiply(fs_vox2ras, v1, v2);
	// v2 contains world coordinates
	MatrixMultiply(fnirt_ras2vox, v2, v3);
	// v3 contains fsl mni voxel coordinates.
	int fslx = V3_X(v3);
	int fsly = V3_Y(v3);
	int fslz = V3_Z(v3);
	GCA_MORPH_NODE *newnode = &(result_gcam->nodes[c][r][s]);

	if ((fslx < 0) || (fslx > (fnirt_atlas->width - 1)) ||
	    (fsly < 0) || (fsly > (fnirt_atlas->height - 1)) ||
	    (fslz < 0) || (fslz > (fnirt_atlas->depth - 1)) ) {
	  // an invalid node
	  newnode->origx = 0;
	  newnode->origy = 0;
	  newnode->origz = 0;

	  newnode->x = 0;
	  newnode->y = 0;
	  newnode->z = 0;
          newnode->invalid = GCAM_POSITION_INVALID;

	} else {
	  GCA_MORPH_NODE *fnirtvcnode = &(fnirt_gcam->nodes[fslx][fsly][fslz]);
	  newnode->origx = fnirtvcnode->origx;
	  newnode->origy = fnirtvcnode->origy;
	  newnode->origz = fnirtvcnode->origz;

	  newnode->x = fnirtvcnode->x;
	  newnode->y = fnirtvcnode->y;
	  newnode->z = fnirtvcnode->z;

	}
      }
    }
  }
  copyVolGeom(&fnirtvolchange_gcam->atlas, &result_gcam->atlas);
  copyVolGeom(&fnirtvolchange_gcam->image, &result_gcam->image);
  
  result_gcam->m_affine = fnirtvolchange_gcam->m_affine;
  result_gcam->det = fnirtvolchange_gcam->det;
  
}
#endif
void writeM3Z(const string& fname, const GCAM *gcam)
// Write an m3z file. Just calls down to GCAMwrite
{
  GCAMwrite(gcam, fname.c_str());
}

int main(int argc, char *argv[])
{
  cout << vcid << endl << endl;

  // Default initialization
  int nargs = handle_version_option(argc, argv, vcid, "$Name:  $");
  if (nargs && argc - nargs == 1)
  {
    exit(0);
  }
  argc -= nargs;
  Progname = argv[0];
  argc--;
  argv++;
  ErrorInit(NULL, NULL, NULL);

  // Parse command line
  if (!parseCommandLine(argc, argv, P))
  {
    //printUsage();
    exit(1);
  }

  GCA_MORPH* gcam = NULL;
  GCA_MORPH* gcamnew = NULL;

  gcam = readM3Z(P.in_warp.c_str());
  if (!gcam)
  {
    ErrorExit(ERROR_BADFILE, "%s: can't read input file %s",
              Progname, P.in_warp.c_str());
  }
  Info(gcam);
  if (P.infoOnly) {
    dumpGCAM(gcam);
    return(0);
  }

  LTA *lta = NULL;
  if (P.affineMat != "") {
    cout << "Loading affine matrix" << endl;
    lta = LTAread(P.affineMat.c_str());
  }
  // Read input transform and convert to RAS2RAS:

  MRI *src = NULL;
  MRI *atlas = NULL;

  if (P.sourceImage != "") {
    src = MRIread(P.sourceImage.c_str());
  }

  if (P.atlasImage != "") {
    atlas=MRIread(P.atlasImage.c_str());
  }

  // allocate a new gcam with the standard atlas size
  GCAM *gcam256 = GCAMalloc(atlas->width, atlas->height, atlas->depth);

  
  resampleGCAMto256(gcam, gcam256);
  // if (P.maskImage != "") {
  //   MRI *mask = MRIread(P.maskImage.c_str());
  //   maskGCAM(gcam256, mask);
  //   MRIfree(&mask);
  // }

  Info(gcam256);
  GCAMfree(&gcam);

  GCA_MORPH* gcam3 = NULL;
  cout << "Downsampling" << endl;
  gcam3 = GCAMdownsample2(gcam256);
  GCAMfree(&gcam256);
  
  if (atlas != NULL && src != NULL) {
    cout << "GCAMchangeVolGeom" << endl;
    /* This is how mri_concatenate_warp changes the source and target */
    gcamnew=GCAMchangeVolGeom(gcam3, src, atlas);
    cout << "GCAMchangeVolGeom - Done" << endl;
  } else {
    ErrorExit(ERROR_BADFILE, "Can't read either %s or %s",
	      P.sourceImage.c_str(), P.atlasImage.c_str());
  }

  if (lta != NULL) {
    if ( src == NULL || atlas == NULL ) {
       cout << "Missing src or destination image, not remapping transform" << endl;
    } else {
      LTAchangeType(lta, LINEAR_RAS_TO_RAS);
      LTAmodifySrcDstGeom(lta, src, atlas);
      LTAchangeType(lta, LINEAR_VOX_TO_VOX);
    }
    //LTAwrite(lta, "mod.lta");
    gcamnew->m_affine = MatrixCopy(lta->xforms[0].m_L, NULL);
    gcamnew->det = MatrixDeterminant(gcamnew->m_affine);
  }

  if (P.maskImage != "") {
    MRI *mask = MRIread(P.maskImage.c_str());
    maskGCAM(gcamnew, mask);
    MRIfree(&mask);
  }


  MRIfree(&src);
  MRIfree(&atlas);

  Info(gcamnew);

  //GCA *gca = GCAread(P.atlasImage.c_str());
 
  //GCAMInitStuff(gcam2, gca, P.mriImage); 
  cout << "Writing" << endl;
  writeM3Z(P.out_warp.c_str(), gcamnew);
  //writeM2Z(P.out_warp.c_str(), gcam);


  GCAMfree(&gcamnew);

  printf("%s successful.\n", Progname);
  return (0);
}

#include "rjb_m3z_stuff.help.xml.h"
static void printUsage(void)
{
  outputHelpXml(rjb_m3z_stuff_help_xml, rjb_m3z_stuff_help_xml_len);
}

/*!
 \fn int parseNextCommand(int argc, char **argv)
 \brief Parses the command-line for next command
 \param   argc  number of command line arguments
 \param   argv  pointer to a character pointer
 \param      P  reference to parameters
 \returns       number of used arguments for this command
 */
static int parseNextCommand(int argc, char *argv[], Parameters & P)
{
  bool have_input = false;
  //bool have_output = false;

  int nargs = 0;
  char *option;

  option = argv[0] + 1;                     // remove '-'
  if (option[0] == '-')
  {
    option = option + 1;  // remove second '-'
  }
  StrUpper(option);

  if (!strcmp(option, "INM3Z") )
  {
    if (have_input) {
      cerr << endl << endl << "ERROR: Only one input warp can be specified"
           << endl << endl;
      printUsage();
      exit(1);
    }
    have_input = true;

    P.in_warp = string(argv[1]);
    P.in_type = filetypes::M3Z;
    nargs = 1;
    cout << "--inm3z: " << P.in_warp << " input M3Z warp." << endl;
  }
  else if (!strcmp(option, "OUTM3Z") )
  {

    P.out_warp = string(argv[1]);
    P.out_type = filetypes::M3Z;
    nargs = 1;
    cout << "--outm3z: " << P.out_warp << " output M3Z." << endl;
  } 
  else if (!strcmp(option, "AFFINEMATRIX") )
  {

    P.affineMat = string(argv[1]);
    nargs = 1;
    cout << "--affineMatrix: " << P.affineMat << " ." << endl;
  } 
  else if (!strcmp(option, "SOURCEIMAGE") )
  {
    P.sourceImage = string(argv[1]);
    nargs=1;
    cout << "--sourceImage: " << P.sourceImage << " ." << endl;
  }
 else if (!strcmp(option, "ATLASIMAGE") )  
  {
    P.atlasImage = string(argv[1]);
    nargs=1;
    cout << "--atlasImage: " << P.atlasImage << " ." << endl;
  }
 else if (!strcmp(option, "MASKIMAGE") ) 
  {
    P.maskImage = string(argv[1]);
    nargs=1;
    cout << "--maskImage: " << P.maskImage << " ." << endl;
  }

  else if (!strcmp(option, "HELP") )
  {
    printUsage();
    exit(1);
  }
  else if (!strcmp(option, "INFO") )
  {
    P.infoOnly = (bool)atoi(argv[1]);
    nargs=1;
    cout << "--info " << P.infoOnly << "." << endl;
  }
  else
  {
    cerr << endl << endl << "ERROR: Option: " << argv[0]
         << " unknown (see --help) !! " << endl << endl;
    exit(1);
  }

  fflush(stdout);

  return (nargs);
}

/*!
 \fn int parseCommandLine(int argc, char **argv)
 \brief Parses the command-line
 \param   argc  number of command line arguments
 \param   argv  pointer to a character pointer
 \param      P  reference to parameters
 \returns       if all necessary parameters were set
 */
static bool parseCommandLine(int argc, char *argv[], Parameters & P)
{
  int nargs;
  int inputargs = argc;
  for (; argc > 0 && ISOPTION(*argv[0]); argc--, argv++)
  {
    nargs = parseNextCommand(argc, argv, P);
    argc -= nargs;
    argv += nargs;
  }

  if (inputargs == 0)
  {
    printUsage();
    exit(1);
  }


  return true;
}
