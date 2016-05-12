#include <string.h>
#include <stdio.h>
#include "fitsio.h"

#include <iostream>
#include <sstream>

#include <sys/time.h>
#include <time.h>
#include <inttypes.h>
#include <fstream>
#include <unistd.h>
#include <getopt.h>    /* for getopt_long; standard getopt is in unistd.h */
#include <vector>
#include <algorithm>
#include <ctime>
#include <climits>
#include <cmath>
#include <iomanip>

#include "globalConstants.h"

using namespace std;

int deleteFile(const char *fileName){
  cout << yellow;
  cout << "Will overwrite: " << fileName << endl << endl;
  cout << normal;
  return unlink(fileName);
}

bool fileExist(const char *fileName){
  ifstream in(fileName,ios::in);
  
  if(in.fail()){
    //cout <<"\nError reading file: " << fileName <<"\nThe file doesn't exist!\n\n";
    in.close();
    return false;
  }
  
  in.close();
  return true;
}

/*========================================================
  ASCII progress bar
==========================================================*/
void showProgress(unsigned int currEvent, unsigned int nEvent) {

  const int nProgWidth=50;

  if ( currEvent != 0 ) {
    for ( int i=0;i<nProgWidth+8;i++)
      cout << "\b";
  }

  double percent = (double) currEvent/ (double) nEvent;
  int nBars = (int) ( percent*nProgWidth );

  cout << " |";
  for ( int i=0;i<nBars-1;i++)
    cout << "=";
  if ( nBars>0 )
    cout << ">";
  for ( int i=nBars;i<nProgWidth;i++)
    cout << " ";
  cout << "| " << setw(3) << (int) (percent*100.) << "%";
  cout << flush;

}

void printCopyHelp(const char *exeName, bool printFullHelp=false){
  
  if(printFullHelp){
    cout << bold;
    cout << endl;
    cout << "This program computes overscan mean and subtracts it line by line.\n";
    cout << "It handles all the available HDUs. The HDUs in the output fit file\n";
    cout << "will be double (64bits) for 64bits images and float (32bis) in all\n";
    cout << "the other cases.\n";
    cout << "The card \"TRIMSEC\" must be present in the header to use this program.\n";
    cout << normal;
  }
  cout << "==========================================================================\n";
  cout << yellow;
  cout << "\nUsage:\n";
  cout << "  "   << exeName << " <input file> -o <output filename> \n\n";
  cout << "\nOptions:\n";
  cout << "  -q for quiet (no screen output)\n";
  cout << "  -f for processing images taken at Fermilab (channels 1 & 12)\n";
  cout << "  -x keep the serial OS in the output image\n";
  cout << "  -y keep the parallel OS in the output image \n";
  cout << "  -s <HDU number> for processing a single HDU \n\n";
  cout << normal;
  cout << blue;
  cout << "For any problems or bugs contact Javier Tiffenberg <javiert@fnal.gov>\n\n";
  cout << normal;
  cout << "==========================================================================\n\n";
}

string bitpix2TypeName(int bitpix){
  
  string typeName;
  switch(bitpix) {
      case BYTE_IMG:
          typeName = "BYTE(8 bits)";
          break;
      case SHORT_IMG:
          typeName = "SHORT(16 bits)";
          break;
      case LONG_IMG:
          typeName = "INT(32 bits)";
          break;
      case FLOAT_IMG:
          typeName = "FLOAT(32 bits)";
          break;
      case DOUBLE_IMG:
          typeName = "DOUBLE(64 bits)";
          break;
      default:
          typeName = "UNKNOWN";
  }
  return typeName;
}


int copyStructure(const string inFile, const char *outF, vector<int> &singleHdu){
    
  fitsfile  *outfptr; /* FITS file pointers defined in fitsio.h */
  fitsfile *infptr;   /* FITS file pointers defined in fitsio.h */
  
  int status = 0;
  int single = 0;
  
  int hdutype, bitpix, naxis = 0, nkeys;
  int nhdu = 0;
  long naxes[9] = {1, 1, 1, 1, 1, 1, 1, 1, 1};
  long totpix = 0;
  //char card[81];
  
  ostringstream fileStructSummary;
  
  const char* inF = inFile.c_str();
  fits_open_file(&infptr, inF, READONLY, &status); /* Open the input file */
  if (status != 0) return(status);
  
  fits_get_num_hdus(infptr, &nhdu, &status);
  
  /* check the extensions to process*/
  for(unsigned int i=0;i<singleHdu.size();++i){
    if(singleHdu[i] > nhdu){
      fits_close_file(infptr,  &status);
      cerr << red << "\nError: the file does not have the required HDU!\n\n" << normal;
      return -1000;
    }
  }
  
  if(singleHdu.size() == 0){
    for(int i=0;i<nhdu;++i){
      singleHdu.push_back(i+1);
    }
  }
  
  const unsigned int nUseHdu = (singleHdu[0]<0)? 1 : singleHdu.size();
  
  if(single)
    fileStructSummary << "The output file will contain " << nUseHdu << " of " << nhdu << " hdus availables in the input files."<< endl;  
  
  fits_create_file(&outfptr, outF, &status);/* Create the output file */
  if (status != 0) return(status);
  
  
  fileStructSummary << bold << "HDU   hdutype  #Axis  #Cols  #Rows   IN_datatype      OUT_datatype\n" << normal;
// HDU  hdutype #Axis #Cols #Rows datatype  
  for (unsigned int eI=0; eI<nUseHdu; ++eI)  /* Main loop through each extension */
  {
    unsigned int n = (singleHdu[0]==-1)? 2 : singleHdu[eI];
    if(singleHdu[0]==-2) n = 1; // Chicago file
    
    /* get image dimensions and total number of pixels in image */
    fits_movabs_hdu(infptr, n, &hdutype, &status);
    for (int i = 0; i < 9; ++i) naxes[i] = 1;
    fits_get_img_param(infptr, 9, &bitpix, &naxis, naxes, &status);
    totpix = naxes[0] * naxes[1] * naxes[2] * naxes[3] * naxes[4] * naxes[5] * naxes[6] * naxes[7] * naxes[8];
    fileStructSummary  << setw (3) << n << "  "  << setw (8) << hdutype << "  " << setw (5) << naxis << "  " << setw (5) << naxes[0] << "  " << setw (5) << naxes[1] << "  " << setw (15) << bitpix2TypeName(bitpix);
    
    //continue;
    if (hdutype != IMAGE_HDU || naxis == 0 || totpix == 0){
      /* just copy tables and null images */
      fits_copy_hdu(infptr, outfptr, 0, &status);
      if (status != 0) return(status);
      fileStructSummary << "  " << setw (8) << magenta << "Not an image HDU" << normal;
    }
    else{
      /* create output image with the same size as the input image*/
      int xMin=0;
      int xMax=naxes[0];
      int yMin=0;
      int yMax=naxes[1];
      fits_get_hdrspace(infptr, &nkeys, NULL, &status);
      char keyName[] = "TRIMSEC";
      char trimsecRecord[1024] = "";
      fits_read_card(infptr, keyName, trimsecRecord, &status);
      if(status==KEY_NO_EXIST){
      	status=0;
      }
      else{
      	string sRec(trimsecRecord);
      	replace(sRec.begin(), sRec.end(), '\'', ' ');
      	replace(sRec.begin(), sRec.end(), '[', ' ');
      	replace(sRec.begin(), sRec.end(), ':', ' ');
      	replace(sRec.begin(), sRec.end(), ',', ' ');
      	replace(sRec.begin(), sRec.end(), ']', ' ');
      	size_t tPos = sRec.find("=");
      	istringstream recISS( sRec.substr(tPos+1) );

      	recISS >> xMin >> xMax >> yMin >> yMax;
        if(gKeepParallelOvsc){
          yMax = naxes[1];
        }
        if(gKeepSerialOvsc == false){
        	naxes[0] = xMax-xMin + 1;
        }
      	naxes[1] = yMax-yMin + 1;
      }
      
      
      if(abs(bitpix) == 64) fits_create_img(outfptr, -64, naxis, naxes, &status);
      else{
        fits_create_img(outfptr, -32, naxis, naxes, &status);
      }
      if (status != 0) return(status);

      /* copy the relevant keywords (not the structural keywords) */
      fits_get_hdrspace(infptr, &nkeys, NULL, &status); 
      for (int i = 1; i <= nkeys; ++i) {
	char card[FLEN_CARD];
        fits_read_record(infptr, i, card, &status);
	if(strncmp (card, "TRIMSEC", 7) == 0) continue;
	if(strncmp (card, "DATASEC", 7) == 0) continue;
        if(strncmp (card, "BZERO", 5) == 0) continue;
        if (fits_get_keyclass(card) > TYP_CMPRS_KEY) fits_write_record(outfptr, card, &status);
      }
      fileStructSummary << "  " << setw (15) << bitpix2TypeName( (abs(bitpix) == 64)? -64:-32 );
    }
    fileStructSummary << endl;
    
  }
  fits_close_file(infptr, &status);
  fits_close_file(outfptr,  &status);

  if(gVerbosity){
    cout << bold << "Files structure summary:\n" << normal;
    cout << fileStructSummary.str();
    cout << green << "Structure copied.\n\n" << normal;
  }
  return status;
}


double mean(const double *v, const int &N){
  
  std::vector<double> temparray(v, v+N);
  std::sort(temparray.begin(), temparray.end());
  
  const int nMin = N/3;
  const int nMax = 2*N/3;
  
  double sum=0;
  for(int i=nMin;i<nMax;++i){
    sum+=temparray[i];
  }
  return sum/(nMax-nMin);
}

void computeOvscMean(fitsfile  *infptr, long *fpixel, long *lpixel, vector<double> &ovrscMean){
  
  int status = 0;
  double nulval = 0.;
  int anynul = 0;
  
  long inc[2]={1,1};
  const int osL  = (lpixel[0]-fpixel[0]+1);
  const int nCol = (lpixel[1]-fpixel[1]+1);
  const long npix =  osL*nCol;
  double* sArray = new double[npix];
  /* Read the images as doubles, regardless of actual datatype. */
  fits_read_subset(infptr, TDOUBLE, fpixel, lpixel, inc, &nulval, sArray, &anynul, &status);

//   cout << osL << " " << nCol << " " << npix << endl;
  
  ovrscMean.resize(nCol);
  for(int l=0;l<nCol;++l){
    ovrscMean[l] = mean(sArray + osL*l, osL ); 
  }
  delete[] sArray;
}


bool isSaturated(const double &pixVal, const int &bitpix, const double &bzero){
  
  float saturationVal = 0;
  
  switch(bitpix) {
      case BYTE_IMG:
          saturationVal = 128+bzero;
          break;
      case SHORT_IMG:
          saturationVal = 32768+bzero;
          break;
      default:
          saturationVal = kSatValue;
  }
  
  if(pixVal>=saturationVal*kSatMargin)
    return true;
  
  return false;
}

void subtractOvscMean(fitsfile  *infptr, long *fpixel, long *lpixel, vector<double> ovrscMean, const int firstColToWrite,const int fullNCol, const int &bitpixIn, const double &bzeroIn, double *outArray, long &nSat){
  
  int status    = 0;
  double nulval = 0.;
  int anynul    = 0;
  
  int nLines     = (lpixel[1]-fpixel[1]+1);    
  const int nCol = (lpixel[0]-fpixel[0]+1); 
  
  
  long inc[2]={1,1};
  const long npix =  nCol*nLines;
  double* sArray = new double[npix];
  /* Read the images as doubles, regardless of actual datatype. */
  fits_read_subset(infptr, TDOUBLE, fpixel, lpixel, inc, &nulval, sArray, &anynul, &status);
  
  // cout << nLines << " " << nCol << " " << npix <<  " " << firstColToWrite <<  " " << fullNCol << endl;
  for(int l=0;l<nLines;++l){
    for(int c=0;c<nCol;++c){
      if( isSaturated(sArray[nCol*l + c], bitpixIn, bzeroIn) ){
        outArray[ fullNCol*l + c + firstColToWrite] = kSatValue;
        ++nSat;
      }
      else{
        outArray[ fullNCol*l + c + firstColToWrite] = sArray[nCol*l + c]-ovrscMean[l];
      }
    }
  }
  
  delete[] sArray;
}


int computeImage(const string inFile, const char *outF, vector<int> &singleHdu){
  int status  = 0;
  int nhduOUT = 0;
  int nhduIN  = 0;
  
  /* Open the output file */
  fitsfile  *outfptr; /* FITS file pointers defined in fitsio.h */
  fits_open_file(&outfptr, outF, READWRITE, &status);
  if (status != 0) return(status);
  fits_get_num_hdus(outfptr, &nhduOUT, &status);
  
  fitsfile  *infptr; /* FITS file pointers defined in fitsio.h */
  const char* inF = inFile.c_str();
  fits_open_file(&infptr, inF, READONLY, &status); /* Open the input file */
  if (status != 0) return(status);
  fits_get_num_hdus(infptr, &nhduIN, &status);
  
  const unsigned int nUseHdu = (singleHdu[0]==-1)? 1 : singleHdu.size();
  const bool isFromChicago = (singleHdu[0]==-2);
  
  ostringstream fileSummary;
  fileSummary << bold << "\nHDU   #SatPixL   #SatPixR    #SatPixTot\n" << normal;
  
  for (unsigned int eI=0; eI<nUseHdu; ++eI)  /* Main loop through each extension */
  {
    
    unsigned int n = (singleHdu[0]==-1)? 2 : singleHdu[eI];
    if(singleHdu[0]==-2) n = 1; // Chicago file
    const int nHDUsToProcess = nUseHdu;
    int nOut = eI+1;
    
    
    /* get output image dimensions and total number of pixels in image */
    int hdutype, bitpix, bytepix, naxis = 0;
    long naxes[9] = {1, 1, 1, 1, 1, 1, 1, 1, 1};
    fits_movabs_hdu(outfptr, nOut, &hdutype, &status);
    for (int i = 0; i < 9; ++i) naxes[i] = 1;
    fits_get_img_param(outfptr, 9, &bitpix, &naxis, naxes, &status);
    long totpix = naxes[0] * naxes[1];
    
    /* Don't try to process data if the hdu is empty */    
    if (hdutype != IMAGE_HDU || naxis == 0 || totpix == 0){
      continue;
    }
    
    bytepix = abs(bitpix) / 8;
    if(bytepix!=4 && bytepix!=8) return -1000;
    
    double* outArray = new double[totpix];
    
    for(int i=0;i<totpix;++i) outArray[i] = 0;
    
    
    /* get input image dimensions and total number of pixels in image */
    int bitpixIn, naxisIn = 0;
    long naxesIn[9] = {1, 1, 1, 1, 1, 1, 1, 1, 1};
    fits_movabs_hdu(infptr, n, &hdutype, &status);
    for (int i = 0; i < 9; ++i) naxesIn[i] = 1;
    fits_get_img_param(infptr, 9, &bitpixIn, &naxisIn, naxesIn, &status);
    //long totpixIn = naxesIn[0] * naxesIn[1];
    double bzeroIn;
    ffgky(infptr, TDOUBLE, "BZERO", &bzeroIn, NULL, &status);
    if (status){
      status = 0;
      bzeroIn = 0.0;
    }
    
    
    if(gVerbosity){
      showProgress((eI)*3+0,3*nHDUsToProcess);
    }
      
    /* Open the input file */
    fits_movabs_hdu(infptr, n, &hdutype, &status);
    int xMin=0;
    int xMax=naxes[0];
    int yMin=0;
    int yMax=naxes[1];
    int nkeys;
    fits_get_hdrspace(infptr, &nkeys, NULL, &status);
    char keyName[] = "TRIMSEC";
    char trimsecRecord[1024] = "";
    fits_read_card(infptr, keyName, trimsecRecord, &status);
    if(status==KEY_NO_EXIST){/* should return an error */
      status=0;
    }
    else{
      string sRec(trimsecRecord);
      replace(sRec.begin(), sRec.end(), '\'', ' ');
      replace(sRec.begin(), sRec.end(), '[', ' ');
      replace(sRec.begin(), sRec.end(), ':', ' ');
      replace(sRec.begin(), sRec.end(), ',', ' ');
      replace(sRec.begin(), sRec.end(), ']', ' ');
      size_t tPos = sRec.find("=");
      istringstream recISS( sRec.substr(tPos+1) );

      recISS >> xMin >> xMax >> yMin >> yMax;
      if(gKeepParallelOvsc){
        yMax = naxes[1];
      }
      if( gKeepSerialOvsc == false){
        naxes[0] = xMax-xMin + 1;
      }
      naxes[1] = yMax-yMin + 1;
    }
    
    long nSatR=0;
    if(!isFromChicago){/* right side */
      if(singleHdu[0]==-1) fits_movabs_hdu(infptr, 2, &hdutype, &status);
      vector<double> ovrscMeanR;
      {
      	long fpixel[2]={1+kPrescan,yMin};
      	long lpixel[2]={xMin-1,yMax};
      	computeOvscMean(infptr, fpixel, lpixel, ovrscMeanR);
      }

      {
      	int xEnd = (xMax + xMin)/2;
      	long fpixel[2]={xMin,yMin};
      	long lpixel[2]={xEnd,yMax};
        if(gKeepSerialOvsc){
          fpixel[0] = 1;   
        }
      	subtractOvscMean(infptr, fpixel, lpixel, ovrscMeanR, 0,naxes[0], bitpixIn, bzeroIn, outArray, nSatR);
      }
    }
    
    if(gVerbosity){
      showProgress((eI)*3+1,3*nHDUsToProcess);
    }
    
    long nSatL=0;
    {/* left side */
      if(singleHdu[0]==-1) fits_movabs_hdu(infptr, nhduIN, &hdutype, &status);
      vector<double> ovrscMeanL;
      {
      	long fpixel[2]={xMax+1,yMin};
      	long lpixel[2]={naxesIn[0]-kPrescan,yMax};
      	computeOvscMean(infptr, fpixel, lpixel, ovrscMeanL);
      }

      {
        long firstColToWrite = naxes[0]/2;
        int xStart = (xMax + xMin)/2 + 1;
        if(isFromChicago){
          xStart = xMin;
          firstColToWrite = 0;
        }
        long fpixel[2]={xStart,yMin};
        long lpixel[2]={xMax,yMax};

        if(gKeepSerialOvsc){
          lpixel[0] = naxes[0];   
        }
      	subtractOvscMean(infptr, fpixel, lpixel, ovrscMeanL, firstColToWrite, naxes[0], bitpixIn, bzeroIn, outArray, nSatL);
      }
    }
    
    if(gVerbosity){
      showProgress((eI)*3+2,3*nHDUsToProcess);
    }

    fits_write_img(outfptr, TDOUBLE, 1, totpix, outArray, &status);
    
    delete[] outArray;
    
    fileSummary  << setw (3) << n << "  "  << setw (9) << nSatL << "  " << setw (9) << nSatR << "  " << setw (12) << nSatL+nSatR << endl;

    if(gVerbosity){
      showProgress((eI)*3+3,3*nHDUsToProcess);
    }
    
  }
  
  /* Close the output file */
  fits_close_file(outfptr,  &status);
  fits_close_file(infptr,  &status);
  if(gVerbosity){
    showProgress(1,1);
    
    cout << endl << fileSummary.str() << endl;
  }
  return status;
}


void checkArch(){
  if(sizeof(float)*CHAR_BIT!=32 || sizeof(double)*CHAR_BIT!=64){
    cout << red;
    cout << "\n ========================================================================================\n";
    cout << "   WARNING: the size of the float and double variables is non-standard in this computer.\n";
    cout << "   The program may malfunction or produce incorrect results\n";
    cout << " ========================================================================================\n";
    cout << normal;
  }
}

int processCommandLineArgs(const int argc, char *argv[], vector<int> &singleHdu, string &inFile, string &outFile){
  
  if(argc == 1) return 1;
  
  bool outFileFlag = false;
  singleHdu.clear();
  int opt=0;
  while ( (opt = getopt(argc, argv, "i:o:s:fcqxyQhH?")) != -1) {
    switch (opt) {
    case 'o':
      if(!outFileFlag){
        outFile = optarg;
        outFileFlag = true;
      }
      else{
        cerr << red << "\nError, can not set more than one output file!\n\n" << normal;
        return 2;
      }
      break;
    case 's':
      singleHdu.push_back(atoi(optarg));
      break;
    case 'f':
      singleHdu.clear();
      singleHdu.push_back(-1);
      cout << yellow << "\nWill process a Fermi file: channels 0 & 11\n" << normal;
      break;
    case 'c':
      singleHdu.clear();
      singleHdu.push_back(-2);
      cout << yellow << "\nWill process a Chicago file: channel 0 & only one side\n" << normal;
      break;
    case 'Q':
    case 'q':
      gVerbosity = 0;
      break;
    case 'y':
      gKeepParallelOvsc = true;
      break;
    case 'x':
      gKeepSerialOvsc = true;
      break;
    case 'h':
    case 'H':
    default: /* '?' */
      return 1;
    }
  }
  
  if(!outFileFlag){
    cerr << red << "\nError: output filename missing.\n" << normal;
    return 2;
  }

  inFile="";
  
  if(argc-optind==0){
    cerr << red << "Error: no input file(s) provided!\n\n" << normal;
    return 1;
  }
  else if(argc-optind>1){
    cerr << red << "Error: more than one input file provided!\n\n" << normal;
    return 1;
  }
  
  inFile=argv[optind];
  if(!fileExist(inFile.c_str())){
    cout << red << "\nError reading input file: " << inFile <<"\nThe file doesn't exist!\n\n" << normal;
    return 1;
  }
  
  return 0;
}

int main(int argc, char *argv[])
{
  
  checkArch(); //Check the size of the double and float variables.
  
  time_t start,end;
  double dif;
  time (&start);
  
  string outFile;
  string inFile;
  vector<int> singleHdu;
  
  int returnCode = processCommandLineArgs( argc, argv, singleHdu, inFile, outFile);
  if(returnCode!=0){
    if(returnCode == 1) printCopyHelp(argv[0],true);
    if(returnCode == 2) printCopyHelp(argv[0]);
    return returnCode;
  }
  
  if(gVerbosity){
    cout << bold << "\nWill read the following file:\n" << normal;
    cout << "\t" << inFile << endl;
    cout << bold << "\nThe output will be saved in the file:\n\t" << normal << outFile << endl;
  }
  
  
  /* Overwrite the output file if it already exist */
  if(fileExist(outFile.c_str())){
    cout << yellow << "\nThe output file exist. " << normal;
    deleteFile(outFile.c_str());
  }
  
  /* Do the actual processing */
  int status = copyStructure( inFile,  outFile.c_str(), singleHdu);
  if (status != 0){ 
    fits_report_error(stderr, status);
    return status;
  }
  
  
  status = computeImage( inFile,  outFile.c_str(), singleHdu);
  if (status != 0){ 
    fits_report_error(stderr, status);
    return status;
  }
  
  /* Report */
  time (&end);
  dif = difftime (end,start);
  if(gVerbosity) cout << green << "All done!\n" << bold << "-> It took me " << dif << " seconds to do it!\n\n" << normal;

  return status;
}
