/* **************************************************************
datph.CPP - object to read and write the structure of pH
               data from/to files
************************************************************** */
using namespace std;
#if !defined(DATPH_H)
  #include "datph.hpp"
#endif

/* *********************************************************** */

Phdata::Phdata(void)
{

  elvend = 1;
  lagpos = -99;
  curpos = 0;

};

/* **************************************************************
                    Public Functions
************************************************************** */

/* *************************************************************
************************************************************* */

int Phdata::get(ifstream& infile)
{

  lagpos = infile.tellg();

  infile >> col >> row;
  infile >> varname;
  infile >> carea;
  infile >> pH;
  infile >> contnent;

  infile.seekg(0, ios::cur);
  curpos = infile.tellg();

  if (curpos < (lagpos + 10)) { elvend = -1;}

  return elvend;

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

int Phdata::getdel(FILE* infile)
{

  elvend = fscanf(infile,"%f,%f, %s ,%d,%lf, %s",
                  &col,&row,varname,&carea,&pH,contnent);
                 
  return elvend;

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

void Phdata::out(ofstream& ofile, float col, float row, char varname[9],
                   int carea, double pH, char contnent[9])
{

  ofile.setf(ios::fixed,ios::floatfield);
  ofile.setf(ios::showpoint);
  ofile.precision(1);

  ofile << col << ' ' << row << ' ';
  ofile << varname << ' ';
  ofile << setprecision(0) << carea << ' ';
  ofile << setprecision(1) << pH << ' ';
  ofile << contnent;
  ofile << endl;

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

void Phdata::outdel(ofstream& ofile, float col, float row, char varname[9],
                      int carea, double pH, char contnent[9])
{

  ofile.setf(ios::fixed,ios::floatfield);
  ofile.setf(ios::showpoint);
  ofile.precision(1);

  ofile << col << "," << row << ", ";
  ofile << varname << " ,";
  ofile << setprecision(0) << carea << ",";
  ofile << setprecision(1) << pH << ", ";
  ofile << contnent;
  ofile << endl;

};

