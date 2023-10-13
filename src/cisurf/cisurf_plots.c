
#include "cisurf.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>




void plotBaseMeshVTK( baseMesh *mesh1, char *filename ) {
  //
  // This routine creates an ASCII vtk file which can plot mesh1.
  //
  
  long nelems, nverts;
  nelems = mesh1->nelems;
  nverts = mesh1->nverts;
  
  baseElement *elem1;
  FILE *fptr;

  fptr = fopen( filename, "w" );
  fprintf(fptr, "# vtk DataFile Version 3.0\n");
  fprintf(fptr, "MESHNAME: %s\n", mesh1->name);
  fprintf(fptr, "ASCII\n");
  fprintf(fptr, "DATASET UNSTRUCTURED_GRID\n");

  fprintf(fptr, "POINTS %ld float\n", nverts);

  long i;
  for (i=0; i<nverts; i++) {
    fprintf(fptr, "%e %e %e\n", mesh1->verts[3*i], mesh1->verts[3*i+1],
	    mesh1->verts[3*i+2]);
  }

  // compute total number of points needed across all elements,
  // allowing for various kinds of elements which have various numbers
  // of vertices defining them
  long ntot;
  ntot = 0;
  for (i=0; i<nelems; i++) {
    ntot = ntot + 1 + mesh1->elements[i].nv;
  }

  fprintf(fptr, "CELLS %ld %ld\n", nelems, ntot);

  long j;
  for (i=0; i<nelems; i++) {
    fprintf(fptr, "%ld", mesh1->elements[i].nv);
    for (j=0; j<mesh1->elements[i].nv; j++) {
      fprintf(fptr, " %ld", mesh1->elements[i].ivs[j]-1);
    }
    fprintf(fptr, "\n");
  }

  // print the cell types
  fprintf(fptr, "CELL_TYPES %ld\n", nelems);
  for (i=0; i<nelems; i++) {
    fprintf(fptr, "22\n");
  }

  // plot the z value of each node so we can color the thing
  fprintf(fptr, "\n");
  fprintf(fptr, "POINT_DATA %ld\n", nverts);
  fprintf(fptr, "SCALARS z-value float 1\n");
  fprintf(fptr, "LOOKUP_TABLE default\n");

  for (i=0; i<nverts; i++) {
    fprintf(fptr, "%e\n", mesh1->verts[3*i+2]);
  }

  /* elseif(ifflat.eq.1) then */
  /*   do i=1,ntri */
  /*     write(iunit1,'(a,i9,i9,i9)') "3 ", Geometry1%Tri(1,i)-1, & */
  /*      Geometry1%Tri(2,i)-1,Geometry1%Tri(3,i)-1 */
  /*   enddo */
  /*   write(iunit1,'(a,i9)') "CELL_TYPES ", ntri */
  /*   do ipatch = 1,ntri */
  /*     write(iunit1,'(a)') "5" */
  /*   end do */
  /* endif */


  fclose(fptr);

  

  return; 
}
