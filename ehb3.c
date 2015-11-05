/*************************************************************************

   Program:    ehb3
   File:       ehb3.c
   
   Version:    V1.0
   Date:       07.02.06
   Function:   Calculate total energy between two H-bonded residues
               specified on the command line
   
   Copyright:  (c) UCL / Dr. Andrew C. R. Martin 2006
   Author:     Dr. Andrew C. R. Martin
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
   V1.0  05.02.03   Original              By: ACRM
   V1.1  07.03.03   Added control file    By: ALC
   V1.2  23.09.05   Fixed various bugs and takes command line parameters
                    for potential type    By: ACRM

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
#define CUTSQ    3.5

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
BOOL gRelax = FALSE;
BOOL gHBOnly = FALSE;

/************************************************************************/
/* Prototypes
*/
BOOL ParseCmdLine(int argc, char **argv, char *PDBFile,
                  char *resspec1, char *resspec2);
BOOL CalcEnergy(char *PDBFile, HBONDS *hbonds, int i);
void Usage(void);
int main(int argc, char **argv);
void FixHydrogenAtomNames(PDB *pdb);
PDB *CopyAndFixResidue(PDB *pdb, char chain);
REAL ParseECalcOutput(char *EnergyFile);
PDB *FixResidue(PDB *res);
PDB *CopyResidue(PDB *pdb, char chain);
BOOL CreateHB(HBONDS *HBonds, char *resspec1, char *resspec2);


/************************************************************************/
/*>int main(int argc, char **argv)
   -------------------------------
   Main program

   06.02.03 Original   By: ACRM
*/
int main(int argc, char **argv)
{
   char    PDBFile[MAXBUFF],
           resspec1[MAXBUFF],
           resspec2[MAXBUFF];
   HBONDS  HBonds[1];
   
   if(ParseCmdLine(argc, argv, PDBFile, resspec1, resspec2))
   {
      if(CreateHB(HBonds, resspec1, resspec2))
      {
         if(!strncmp(HBonds[0].type, "SS", 2) &&
            strncmp(HBonds[0].AtomD, "OXT", 3) &&
            strncmp(HBonds[0].AtomA, "OXT", 3))
         {
            if(!CalcEnergy(PDBFile, HBonds, 0))
               return(1);
         }
         else
         {
            printf("Error: Can't calc energy involving non-SS or OXT\n");
            return(1);
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
   23.09.05 Updated for V1.2
*/
void Usage(void)
{
   fprintf(stderr,"\nehb3 V3 (c) 2003-6, Dr. Andrew C.R. Martin, UCL, The \
University of Reading\n");
   fprintf(stderr, "and Dr. Alison L. Cuff\n");
   
   fprintf(stderr,"\nUsage: ehb3 [-r][-o] pdhfile resspec1 resspec2\n");
   fprintf(stderr,"\n       -r  Use the RELAX option in ecalc\n");
   fprintf(stderr,"       -o  Calculate the hbond energy only (overrides -r)\n");

   fprintf(stderr,"\n       pdhfile    - PDB file with hydrogens\n");
   fprintf(stderr,"       resspec - residue and atom specifier in the form \
[c]nnn[i].atom\n");

   fprintf(stderr,"\nehb3 calculates the total energy for a pair of amino \
acids identified\n");
   fprintf(stderr,"as being in a sidechain-sidechain hydrogen bond\n\n");
}

/************************************************************************/
/*>BOOL ParseCmdLine(int argc, char **argv, char *PDBFile, 
                     char *resspec1, char *resspec2)
   -------------------------------------------------------
   Parse the command line

   06.02.03 Original   By: ACRM
*/
BOOL ParseCmdLine(int argc, char **argv, char *PDBFile,
                  char *resspec1, char *resspec2)
{
   argc--;
   argv++;
   
   while(argc && argv[0][0] == '-')
   {
      switch(argv[0][1])
      {
      case 'r':
         gRelax = TRUE;
         break;
      case 'o':
         gHBOnly = TRUE;
         break;
      case 'h':
         return(FALSE);
      default:
         return(FALSE);
      }
      
      argc--;
      argv++;
   }
   
   if(argc != 3)
      return(FALSE);
   
   strcpy(PDBFile,argv[0]);
   strcpy(resspec1,argv[1]);
   strcpy(resspec2,argv[2]);
   
   return(TRUE);
}


/************************************************************************/
/*>BOOL CalcEnergy(char *PDBFile, HBONDS *hbonds,  int i)
   -----------------------------------------------------
   Given a PDB file with hydrogens (in Charmm format) and a list of
   HBonds, calculate the energy for one HBond

   06.02.03 Original   By: ACRM
*/
BOOL CalcEnergy(char *PDBFile, HBONDS *hbonds, int i)
{
   FILE       *fp;
   char       chainA[8],  chainD[8],
              insertA[8], insertD[8],
              cmdline[MAXBUFF],
              PDBFilename[MAXBUFF],
              EnergyFile[MAXBUFF],
              CONTROLfile[MAXBUFF];
   
   int        resnumA, resnumD,
              natoms;
   REAL       energy;
   static PDB *pdb = NULL;
   PDB        *donor,
              *acceptor,
              *p,
              *acceptor_c,
              *donor_c,
              *acceptor_n,
              *donor_n;
   
   sprintf(CONTROLfile, "control.dat.%d", (int)getpid());
   
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
/*         FixHydrogenAtomNames(pdb); */
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
   ParseResSpec(hbonds[i].ResID_A, chainA, &resnumA, insertA);
   ParseResSpec(hbonds[i].ResID_D, chainD, &resnumD, insertD);

   /* Find these residues                                               */
   if((donor = FindResidue(pdb, chainD[0], resnumD, insertD[0]))==NULL)
   {
      fprintf(stderr,"Donor residue %s not found\n", 
              hbonds[i].ResID_D);
      return(FALSE);
   }
   
   if((acceptor = FindResidue(pdb, chainA[0], resnumA, insertA[0]))==NULL)
   {
      fprintf(stderr,"Acceptor residue %s not found\n", 
              hbonds[i].ResID_A);
      return(FALSE);
   }

   /* 22.09.05 If the two residues are bonded, join them into one       */
   acceptor_c = FindAtomInRes(acceptor, "C   ");
   donor_c    = FindAtomInRes(donor,    "C   ");
   acceptor_n = FindAtomInRes(acceptor, "N   ");
   donor_n    = FindAtomInRes(donor,    "N   ");
   if(DISTSQ(acceptor_c, donor_n) < CUTSQ)
   {
      if((acceptor = CopyResidue(acceptor, 'X'))==NULL)
         return(FALSE);
      if((donor = CopyResidue(donor, 'X'))==NULL)
         return(FALSE);
      
      p = acceptor;
      LAST(p);
      p->next = donor;
      acceptor = FixResidue(acceptor);
      donor = NULL;
   }
   else if(DISTSQ(donor_c, acceptor_n) < CUTSQ)
   {
      if((acceptor = CopyResidue(acceptor, 'X'))==NULL)
         return(FALSE);
      if((donor = CopyResidue(donor, 'X'))==NULL)
         return(FALSE);
      
      p = donor;
      LAST(p);
      p->next = acceptor;
      donor = FixResidue(donor);
      acceptor = NULL;
   }
   else
   {
      /* Copy just these residues, add NTER and CTER and fix chain names*/
      if((acceptor = CopyAndFixResidue(acceptor, 'A'))==NULL)
         return(FALSE);
      if((donor    = CopyAndFixResidue(donor, 'D'))==NULL)
         return(FALSE);
   }
   
   /* Now write temporary PDB file containing just donor and acceptor 
      residues 
   */
          
   sprintf(PDBFilename, "%d.pdh", (int)getpid());
   sprintf(EnergyFile,  "%d.ec",  (int)getpid());

   /* Write PDBFilename and other options  to control file */
   if((fp=fopen(CONTROLfile, "w"))!=NULL)
   {    
      fprintf(fp, "PDBFILE %s\n", PDBFilename);
      fprintf(fp, "IGNTER\n");
      if(gHBOnly)
      {
         fprintf(fp, "POTENTIAL\n");
         fprintf(fp, "HBONDS\n");
         fprintf(fp, "END\n");
      }
      else if(gRelax)
      {
         fprintf(fp, "RELAX\n");
      }
      fclose(fp);
   }
   else
   {
      fprintf(stderr,"Can't write control file: %s\n", CONTROLfile);
      return(FALSE);
   }
   

   /* Write donor/acceptor residues to temporary PDB file */

   if((fp=fopen(PDBFilename, "w"))!=NULL)
   {
      for(p=donor; p!=NULL; NEXT(p))
         WritePDBRecordAtnam(fp, p);
      for(p=acceptor; p!=NULL; NEXT(p))
         WritePDBRecordAtnam(fp, p);
      fclose(fp);

      /* Call the ecalc program to calculate the energy  */
 
      sprintf(cmdline, "ecalc %s > %s\n", CONTROLfile, EnergyFile);
       
      system(cmdline); 

      /* Extract the energy from the output of ecalc                    */
      if((energy = ParseECalcOutput(EnergyFile))==(REAL)-99999.999)
         return(FALSE);

      /* Print the energy                                               */
      fprintf(stdout, "%.6f\n", energy);
      
      unlink(PDBFilename);
      unlink(EnergyFile);
   }
   else
   {
      fprintf(stderr,"Can't write temporary PDB file\n");
      return(FALSE);
   }

   /* Free the copies of the residues                                   */
   FREELIST(donor, PDB);
   FREELIST(acceptor, PDB);

   /* Delete the control file                                           */
   unlink(CONTROLfile);
   
   return(TRUE);
}

/************************************************************************/
/*>PDB *CopyResidue(PDB *pdb, char chain)
   --------------------------------------
   Creates a copy of a residue
   Returns a new linked list

   22.09.05 Original   By: ACRM
*/
PDB *CopyResidue(PDB *pdb, char chain)
{
   PDB *p, *r=NULL, 
       *end, 
       *res = NULL;
   
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

   return(res);
}

/************************************************************************/
/*>PDB *FixResidue(PDB *res)
   -------------------------
   Adds nter and cter residues and atoms
   (Should test the called routines!)

   22.09.05 Original   By: ACRM
*/
PDB *FixResidue(PDB *res)
{
   /* Now add the NTER residue                                          */
   AddNTerHs(&res, TRUE); 
   /* Now add the CTER residue                                          */
   FixCterPDB(res, 2); 

   return(res);
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
   PDB *res = NULL;
   
   res = CopyResidue(pdb, chain);
   res = FixResidue(res);

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
         p->atnam_raw[4] = '\0';
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

/************************************************************************/
/*>int ReadHBonds(char *filename, HBONDS *HBonds)
   ----------------------------------------------
   Reads the HBond list from HBPlus output

   Lifted from ehb.c

   04.01.95 Original    By: ACRM  (from ehb.c)
   05.02.03 Modified to store residue ID and name as well
*/
BOOL CreateHB(HBONDS *HBonds, char *resspec1, char *resspec2)
{
   char *stop;
   
   UPPER(resspec1);
   UPPER(resspec2);
   
   if((stop = strchr(resspec1, '.'))==NULL)
      return(FALSE);
   *stop = '\0';
   strncpy(HBonds[0].ResID_D, resspec1, 8);
   strncpy(HBonds[0].AtomD, stop+1, 8);
   
   if((stop = strchr(resspec2, '.'))==NULL)
      return(FALSE);
   *stop = '\0';
   strncpy(HBonds[0].ResID_A, resspec2, 8);
   strncpy(HBonds[0].AtomA, stop+1, 8);

   strcpy(HBonds[0].type, "SS");

   return(TRUE);
}

