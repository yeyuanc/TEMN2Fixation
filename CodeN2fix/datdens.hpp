/* **************************************************************
DATDENS.HPP - object to read and write the structure of density
               data from/to files
************************************************************** */
using namespace std;
#if !defined(DATDENS_H)
#define DATDENS_H

class Densdata {

  public:

     Densdata(void);

/* **************************************************************
                      Public Functions
************************************************************** */

// read data structure.
     int get(ifstream& infile);
     int getdel(FILE* infile);
//write data structure.
     void out(ofstream& ofile, float col, float row, char varname[9], int carea, double density, char contnent[9]);
     void outdel(ofstream& ofile, float col, float row, char varname[9], int carea, double density, char contnent[9]);


/* **************************************************************
                     Public Variables
************************************************************** */

     float col;          // column or longitude of grid cell (degrees)
     float row;          // row or latitude of grid cell (degrees)
     char varname[9];    // "ELEV"
     int carea;          // area covered by grid cell (sq. km)
     double density;        // density of grid cell
     char contnent[9];   // name of continent containing grid cell


  private:

/* **************************************************************
                      Private Variables
************************************************************** */

     int elvend;
     long curpos;
     long lagpos;

};

#endif

