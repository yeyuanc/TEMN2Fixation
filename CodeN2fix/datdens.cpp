/* **************************************************************
datph.CPP - object to read and write the structure of pH
               data from/to files
************************************************************** */
using namespace std;
#if !defined(DATDENS_H)
  #include "datdens.hpp"
#endif

/* *********************************************************** */

Densdata::Densdata(void)
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

int Densdata::get(ifstream& infile)
{

  lagpos = infile.tellg();

  infile >> col >> row;
  infile >> varname;
  infile >> carea;
  infile >> density;
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

int Densdata::getdel(FILE* infile)
{

  elvend = fscanf(infile,"%f,%f, %s ,%d,%lf, %s",
                  &col,&row,varname,&carea,&density,contnent);
                 
  return elvend;

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

void Densdata::out(ofstream& ofile, float col, float row, char varname[9],
                   int carea, double density, char contnent[9])
{

  ofile.setf(ios::fixed,ios::floatfield);
  ofile.setf(ios::showpoint);
  ofile.precision(1);

  ofile << col << ' ' << row << ' ';
  ofile << varname << ' ';
  ofile << setprecision(0) << carea << ' ';
  ofile << setprecision(1) << density << ' ';
  ofile << contnent;
  ofile << endl;

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

void Densdata::outdel(ofstream& ofile, float col, float row, char varname[9],
                      int carea, double density, char contnent[9])
{

  ofile.setf(ios::fixed,ios::floatfield);
  ofile.setf(ios::showpoint);
  ofile.precision(1);

  ofile << col << "," << row << ", ";
  ofile << varname << " ,";
  ofile << setprecision(0) << carea << ",";
  ofile << setprecision(1) << density << ", ";
  ofile << contnent;
  ofile << endl;

};

