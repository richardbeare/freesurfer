/**
 * @file  rjb_m3z_stuff.cpp
 * @brief A program to explore warp files
 *
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
  filetypes::FileType in_type;
  filetypes::FileType out_type;
};

static struct Parameters P =
  { "", "", "", filetypes::UNKNOWN, filetypes::UNKNOWN};

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
  int row;

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

  for (row=1; row < 5; ++row)
  {
    printf("              %8.4f %8.4f %8.4f %10.4f\n",
           mat->rptr[row][1], mat->rptr[row][2],
           mat->rptr[row][3], mat->rptr[row][4]);
  }

  return (NO_ERROR);
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

void Info(const GCAM* gcam)
{
  cout << "Width, height, depth " << gcam->width << " " << gcam->height << " " << gcam->depth << endl;
  cout << "Atlas VG" << endl;
  VGInfo(&(gcam->atlas));
  cout << endl;
  cout << "Image VG" << endl;
  VGInfo(&(gcam->image));
  PrettyMatrixPrint((gcam->m_affine));
  
  cout << "Type : " << gcam->type << endl;
  cout << "status : " << gcam->status << endl;

  
}

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

  LTA *lta = NULL;
  if (P.affineMat != "") {
    cout << "Loading affine matrix" << endl;
    lta = LTAread(P.affineMat.c_str());
  }
  // Read input transform and convert to RAS2RAS:
  GCA_MORPH* gcam = NULL;
  gcam = readM3Z(P.in_warp.c_str());

  if (!gcam)
  {
    ErrorExit(ERROR_BADFILE, "%s: can't read input file %s",
              Progname, P.in_warp.c_str());
  }

  if (lta != NULL) {
    gcam->m_affine = lta->xforms->m_L;
  }

  if (gcam->m_affine == NULL) {
    cout << "Initializing identity matrix" << endl;
    gcam->m_affine = MatrixIdentity(4, NULL);
  }

  Info(gcam);

  cout << "Downsampling" << endl;

  GCA_MORPH* gcam2 = NULL;
  gcam2 = GCAMdownsample2(gcam);

  cout << "Writing" << endl;
  writeM3Z(P.out_warp.c_str(), gcam2);
  //writeM3Z(P.out_warp.c_str(), gcam);


  GCAMfree(&gcam);
  GCAMfree(&gcam2);

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
  else if (!strcmp(option, "HELP") )
  {
    printUsage();
    exit(1);
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
