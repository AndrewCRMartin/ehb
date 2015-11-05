/*************************************************************************

   Program:    ehb2
   File:       ehb2.c
   
   Version:    V1.0
   Date:       05.02.03
   Function:   Calculate total energy between two H-bonded residues
               identified by HBPlus
   
   Copyright:  (c) University of Reading / Dr. Andrew C. R. Martin 2003
   Author:     Dr. Andrew C. R. Martin
   Address:    School of Animal and Microbial Sciences,
               The University of Reading,
               Whiteknights,
               P.O. Box 228,
               Reading RG6 6AJ.
               England.
   Phone:      +44 (0)118 987 5123 Extn. 7022
   Fax:        +44 (0)118 931 0180
   EMail:      andrew@bioinf.org.uk
               
**************************************************************************

   This program is not in the public domain, but it may be copied
   according to the conditions laid out in the accompanying file
   COPYING.DOC

   The code may be modified as required, but any modifications must be
   documented so that the person responsible can be identified. If someone
   else breaks this code, I don't want to be blamed for code that does not
   work! 

   The code may not be sold commercially or included as part of a 
   commercial product except as described in the file COPYING.DOC.

**************************************************************************

   Description:
   ============

**************************************************************************

   Usage:
   ======

**************************************************************************

   Revision History:
   =================

*************************************************************************/
/* Includes
*/
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include "bioplib/macros.h"
#include "bioplib/pdb.h"
#include "bioplib/fsscanf.h"
#include "bioplib/MathUtil.h"

/************************************************************************/
/* Defines and macros
*/
#define MAXHBOND 10000
#define MAXBUFF  256
#define NSKIP    8  /* Number of header lines at start of HBPlus output */

typedef struct
{
   char AtomD[8],
        AtomA[8],
        AtomH[8],
        ResID_D[8],
        ResID_A[8],
        Resnam_D[8],
        Resnam_A[8],
        type[8];
   REAL DistDA,
        AngDHA,
        DistHA,
        AngHAAA,
        AngDAAA;
}  HBONDS;

/************************************************************************/
/* Globals
*/

/************************************************************************/
/* Prototypes
*/
BOOL ParseCmdLine(int argc, char **argv, char *PDBFile, char *HBPlusFile,
                  BOOL *relax, BOOL *all);
BOOL CalcEnergy(char *PDBFile, HBONDS *hbonds, int i, BOOL relax);
void Usage(void);
int main(int argc, char **argv);
int ReadHBonds(char *filename, HBONDS *HBonds);
void FixHydrogenAtomNames(PDB *pdb);
PDB *CopyAndFixResidue(PDB *pdb, char chain);
REAL ParseECalcOutput(char *EnergyFile);


/************************************************************************/
/*>int main(int argc, char **argv)
   -------------------------------
   Main program

   06.02.03 Original   By: ACRM
*/
int main(int argc, char **argv)
{
   char    PDBFile[MAXBUFF],
           HBPlusFile[MAXBUFF];
   HBONDS  HBonds[MAXHBOND];
   int     NHBonds, i;
   BOOL    relax,
           all;
   
   
   if(ParseCmdLine(argc, argv, PDBFile, HBPlusFile, &relax, &all))
   {
      if((NHBonds = ReadHBonds(HBPlusFile, HBonds))==0)
         return(1);
      
      for(i=0; i<NHBonds; i++)
      {
         if(all)
         {
            if(!CalcEnergy(PDBFile, HBonds, i, relax))
               return(1);
         }
         else
         {
            /* If this is a sidechain-sidechain HBonds and doesn't involve
               the OXT atoms (a bug in HBPlus...), then calculate the 
               energy
            */
            if(!strncmp(HBonds[i].type, "SS", 2) &&
               strncmp(HBonds[i].AtomD, "OXT", 3) &&
               strncmp(HBonds[i].AtomA, "OXT", 3))
            {
               if(!CalcEnergy(PDBFile, HBonds, i, relax))
                  return(1);
            }
         }
      }
   }
   else
   {
      Usage();
   }

   return(0);
}

/************************************************************************/
/*>void Usage(void)
   ----------------
   Print Usage message 

   06.02.03 Original   By: ACRM
*/
void Usage(void)
{
   fprintf(stderr,"\nehb2 V1.0 (c) 2003, Dr. Andrew C.R. Martin, The \
University of Reading\n");

   fprintf(stderr,"\nUsage: ehb2 pdhfile hbplusfile\n");

   fprintf(stderr,"\n       pdhfile    - PDB file with hydrogens output \
from HBPlus (xxxx.h)\n");
   fprintf(stderr,"       hbplusfile - the main results file from HBPlus \
(xxxx.hb2)\n");

   fprintf(stderr,"\nehb2 calculates the total energy for a pair of amino \
acids identified\n");
   fprintf(stderr,"as being in a sidechain-sidechain hydrogen bond by \
HBPlus. HBPlus\n");
   fprintf(stderr,"should be run with the -o flag in order to generate a \
PDB file \n");
   fprintf(stderr,"containing hydrogens. This is used as input to ehb2 \
together with the\n");
   fprintf(stderr,"main HBPlus output file.\n\n");
}

/************************************************************************/
/*>BOOL ParseCmdLine(int argc, char **argv, char *PDBFile, 
                     char *HBPlusFile, BOOL *relax, BOOL *all)
   -----------------------------------------------------------
   Parse the command line

   06.02.03 Original   By: ACRM
   07.02.03 Added relax parameter (-r) and
            Added all parameter (-a)
*/
BOOL ParseCmdLine(int argc, char **argv, char *PDBFile, char *HBPlusFile,
                  BOOL *relax, BOOL *all)
{
   *relax = FALSE;
   *all   = FALSE;

   argc--;
   argv++;
   
   while(argc && argv[0][0] == '-')
   {
      switch(argv[0][1])
      {
      case 'a':
         *all = TRUE;
         break;
      case 'r':
         *relax = TRUE;
         break;
      default:
         return(FALSE);
      }
      
      argc--;
      argv++;
   }
   
   if(argc != 2)
      return(FALSE);
   
   strcpy(PDBFile,argv[0]);
   strcpy(HBPlusFile,argv[1]);

   return(TRUE);
}

/************************************************************************/
/*>int ReadHBonds(char *filename, HBONDS *HBonds)
   ----------------------------------------------
   Reads the HBond list from HBPlus output

   Lifted from ehb.c

   04.01.95 Original    By: ACRM  (from ehb.c)
   05.02.03 Modified to store residue ID and name as well
   07.02.03 Added check for blank lines
*/
int ReadHBonds(char *filename, HBONDS *HBonds)
{
   FILE *fp     = NULL;
   int  NHBonds = 0,
        i;
   char buffer[MAXBUFF];
   
   /* Open the file for reading                                         */
   if((fp=fopen(filename,"r"))!=NULL)
   {
      /* Skip the first NSKIP lines                                     */
      for(i=0; i<NSKIP; i++)
         fgets(buffer,MAXBUFF,fp);

      while(fgets(buffer,MAXBUFF,fp))
      {
         fsscanf(buffer,"%6s%3s%1x%3s%1x%6s%3s%1x%3s%5lf%1x%2s%10x%6lf%1x%5lf%6lf%6lf",
                 HBonds[NHBonds].ResID_D,
                 HBonds[NHBonds].Resnam_D,
                 HBonds[NHBonds].AtomD,
                 HBonds[NHBonds].ResID_A,
                 HBonds[NHBonds].Resnam_A,
                 HBonds[NHBonds].AtomA,
                 &(HBonds[NHBonds].DistDA),
                 HBonds[NHBonds].type,
                 &(HBonds[NHBonds].AngDHA),
                 &(HBonds[NHBonds].DistHA),
                 &(HBonds[NHBonds].AngHAAA),
                 &(HBonds[NHBonds].AngDAAA));

         if(strncmp(HBonds[NHBonds].ResID_D, "      ", 6))
         {
            HBonds[NHBonds].AngDHA  *= PI / (REAL)180.0;
            HBonds[NHBonds].AngHAAA *= PI / (REAL)180.0;
            HBonds[NHBonds].AngDAAA *= PI / (REAL)180.0;
            
            if((++NHBonds) >= MAXHBOND)
            {
               fprintf(stderr,"Too many HBonds, Increase MAXHBONDS\n");
               return(0);
            }
         }
      }

      fclose(fp);
   }
   
   return(NHBonds);
}


/************************************************************************/
/*>BOOL CalcEnergy(char *PDBFile, HBONDS *hbonds, int i, BOOL relax)
   -----------------------------------------------------------------
   Given a PDB file with hydrogens (in Charmm format) and a list of
   HBonds, calculate the energy for one HBond

   06.02.03 Original   By: ACRM
   07.02.03 Added relax parameter
*/
BOOL CalcEnergy(char *PDBFile, HBONDS *hbonds, int i, BOOL relax)
{
   FILE       *fp;
   char       chainA,  chainD,
              insertA, insertD,
              cmdline[MAXBUFF],
              PDBFilename[MAXBUFF],
              ECFile[MAXBUFF],
              EnergyFile[MAXBUFF];
   int        resnumA, resnumD,
              natoms;
   REAL       energy;
   static PDB *pdb = NULL;
   PDB        *donor,
              *acceptor,
              *p;
   

   /* If the PDB file hasn't been read in yet, then do so               */
   if(pdb == NULL)
   {
      if((fp=fopen(PDBFile, "r"))!=NULL)
      {
         if((pdb = ReadPDB(fp, &natoms))==NULL)
         {
            fprintf(stderr,"Can't read atoms from PDB file %s\n",
                    PDBFile);
            fclose(fp);
            return(FALSE);
         }
         FixHydrogenAtomNames(pdb);
         fclose(fp);
      }
      else
      {
         fprintf(stderr,"Unable to open PDB file with hydrogens: %s\n",
                 PDBFile);
         return(FALSE);
      }
   }
   
   /* Find the details of the donor and acceptor residues               */
   fsscanf(hbonds[i].ResID_A, "%c%4d%c", &chainA, &resnumA, &insertA);
   fsscanf(hbonds[i].ResID_D, "%c%4d%c", &chainD, &resnumD, &insertD);
   if(chainA == '-')      chainA  = ' ';
   if(chainD == '-')      chainD  = ' ';
   if(insertA == '-')     insertA = ' ';
   if(insertD == '-')     insertD = ' ';

   /* Find these residues                                               */
   if((donor = FindResidue(pdb, chainD, resnumD, insertD))==NULL)
   {
      fprintf(stderr,"Donor residue %s not found\n", 
              hbonds[i].ResID_D);
      return(FALSE);
   }
   
   if((acceptor = FindResidue(pdb, chainA, resnumA, insertA))==NULL)
   {
      fprintf(stderr,"Acceptor residue %s not found\n", 
              hbonds[i].ResID_A);
      return(FALSE);
   }

   /* Copy just these residues, add NTER and CTER and fix chain names   */
   if((acceptor = CopyAndFixResidue(acceptor, 'A'))==NULL)
      return(FALSE);
   if((donor    = CopyAndFixResidue(donor, 'D'))==NULL)
      return(FALSE);
   
   /* Now write temporary PDB file containing these residues            */
   sprintf(PDBFilename, "/tmp/%d.pdh", (int)getpid());
   sprintf(EnergyFile,  "/tmp/%d.ec",  (int)getpid());
   if((fp=fopen(PDBFilename, "w"))!=NULL)
   {
      for(p=donor; p!=NULL; NEXT(p))
         WritePDBRecord(fp, p);
      for(p=acceptor; p!=NULL; NEXT(p))
         WritePDBRecord(fp, p);
      fclose(fp);

      FREELIST(donor, PDB);
      FREELIST(acceptor, PDB);

      /* If we want to relax the structure we need to write a control
         file for ECalc
      */
      if(relax)
      {
         sprintf(ECFile,  "/tmp/%d.ecalc",  (int)getpid());
         if((fp=fopen(ECFile, "w"))!=NULL)
         {
            fprintf(fp,"PDBFILE %s\n", PDBFilename);
            fprintf(fp,"RELAX\n");
            fclose(fp);
         }
         else
         {
            unlink(PDBFilename);
            unlink(EnergyFile);
            fprintf(stderr,"Can't write ECalc control file\n");
            return(FALSE);
         }
      }
      

      /* Call the ecalc program to calculate the energy                 */
      if(relax)
      {
         sprintf(cmdline, "ecalc %s > %s\n", ECFile, EnergyFile);
      }
      else
      {
         sprintf(cmdline, "ecalc -p %s > %s\n", PDBFilename, EnergyFile);
      }
      system(cmdline);

      /* Extract the energy from the output of ecalc                    */
      if((energy = ParseECalcOutput(EnergyFile))==(REAL)-99999.999)
         return(FALSE);

      /* Print the energy                                               */
      fprintf(stdout, "HBond %d Energy: %.6f\n", i+1, energy);
      
      if(relax) unlink(ECFile);
      unlink(PDBFilename);
      unlink(EnergyFile);
   }
   else
   {
      fprintf(stderr,"Can't write temporary PDB file\n");
      return(FALSE);
   }

   return(TRUE);
}

/************************************************************************/
/*>PDB *CopyAndFixResidue(PDB *pdb, char chain)
   --------------------------------------------
   Make a copy of the PDB linked list for a single residue and add
   NTER and CTER residues to that. Fix the chain name for all to that
   given.

   Returns a new PDB linked list

   06.02.03 Original   By: ACRM
*/
PDB *CopyAndFixResidue(PDB *pdb, char chain)
{
   PDB *end,
       *res = NULL,
       *p,
       *r;
   
   /* Create a copy of the residue                                      */
   end = FindNextResidue(pdb);
   for(p=pdb; p!=end; NEXT(p))
   {
      if(res==NULL)
      {
         INIT(res, PDB);
         r = res;
      }
      else
      {
         ALLOCNEXT(r, PDB);
      }
      if(r==NULL)
      {
         fprintf(stderr,"No memory to copy residue");
         return(NULL);
      }
      
      *r = *p;
      r->chain[0] = chain;
      r->chain[1] = '\0';
      r->next = NULL;
   }

   /* Now add the NTER residue                                          */
   AddNTerHs(&res, TRUE); 
   /* Now add the CTER residue                                          */
   FixCterPDB(res, 2); 

   return(res);
}


/************************************************************************/
/*>void FixHydrogenAtomNames(PDB *pdb)
   -----------------------------------
   For hydrogen atoms, copy the name we have created across to the
   raw name that will be printed out

   06.02.03 Original   By: ACRM
*/
void FixHydrogenAtomNames(PDB *pdb)
{
   PDB *p;
   
   for(p=pdb; p!=NULL; NEXT(p))
   {
      if(p->atnam[0] == 'H')
      {
         strcpy(p->atnam_raw, " ");
         strcat(p->atnam_raw, p->atnam);
      }
   }
}


/************************************************************************/
/*>REAL ParseECalcOutput(char *EnergyFile)
   ---------------------------------------
   Grab the energy from the output file generated by ECalc. Returns
   the magic number -99999.999 on failure.

   06.02.03 Original   By: ACRM
*/
REAL ParseECalcOutput(char *EnergyFile)
{
   FILE *fp;
   char buffer[MAXBUFF];
   REAL energy;
   
   if((fp=fopen(EnergyFile, "r"))!=NULL)
   {
      fgets(buffer, MAXBUFF, fp);
      fgets(buffer, MAXBUFF, fp);
      sscanf(buffer,"%*s %*s %*s %*s %lf", &energy);
      return(energy);
   }
   
   return((REAL)-99999.999);
}

