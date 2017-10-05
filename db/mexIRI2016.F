#include "fintrf.h"
C======================================================================
C     mex -g -v -output mexIRI2016 *.f
C======================================================================
C mexIRI2016.mexw64 was compiled by Alexei Skazik (e-mail: aiskazik@sfedu.ru).
C======================================================================
#if 0
C     
C     mexIRI2016.F
C     .F file needs to be preprocessed to generate .for equivalent
C     
#endif
C     
C     mexIRI2016.f
C
C     Gateway function for irisub.f
      
C     This is a MEX-file for MATLAB.
C     Copyright 1984-2011 The MathWorks, Inc.
C     16.04.2011
C     $Revision: 0.0.0.1 $
C======================================================================
C     Gateway routine
      subroutine mexFunction(nlhs, plhs, nrhs, prhs)

C     Declarations
      implicit none

C     mexFunction arguments:
      mwPointer plhs(*), prhs(*)
      integer nlhs, nrhs

C     Function declarations:
      mwPointer mxGetPr
      mwPointer mxCreateNumericMatrix
      integer mxIsNumeric
      mwPointer mxGetM, mxGetN
      integer mxClassIDFromClassName

C     Pointers to input/output mxArrays:
      mwPointer P_OUTF, P_OARR
      real outf(20,1000),oarr(100,1000)
      logical jf(50)

C     Array information:
      mwPointer mrows, ncols
      mwSize size

C     Arguments for computational routine:
      integer M, N, j
      integer JMAG
      real*8 JF_r(50)
      real*8 oarr_r(100)
      real*8 ALATI, ALONG
      integer IYYYY, MMDD, IUT, HEIGHT
      real*8 DHOUR, H_TEC_MAX 
      integer IVAR
      real*8 VBEG, VEND, VSTP, x

C-----------------------------------------------------------------------
C     Check for proper number of arguments. 
      if(nrhs .ne. 15) then
         call mexErrMsgIdAndTxt ('MATLAB:mexIRI2016:nInput',
     + 'mexIRI2016 requires 15 input arguments.')
      elseif(nlhs .gt. 2) then
         call mexErrMsgIdAndTxt ('MATLAB:mexIRI2016:nOutput',
     +'Too many output arguments. mexIRI2016 requires 2 ',
     +'output arguments.')
      endif

C     Validate inputs
C     Check that the input is a number.
      do j = 1, 15
        if(mxIsNumeric(prhs(j)) .eq. 0) then
            call mexErrMsgIdAndTxt ('MATLAB:mexIRI2016:NonNumeric',
     +                              'Input must be a number.')
        endif
      end do

C     Get the size of the input array JF. IT CAN BE 50 X 1 OR 1 X 50.
      mrows = mxGetM(prhs(2))
      ncols = mxGetN(prhs(2))
      size = mrows*ncols
      
      if (size .ne. 50) then
        call mexErrMsgIdAndTxt('mexIRI2016 requires that JF be',
     +' a 50 x 1 vector.')
      endif

C     Create Fortran data from input arguments.
c     Все входные и выходные аргументы типа REAL*8 = DOUBLE
        
        call mxCopyPtrToReal8(mxGetPr(PRHS(1)),x,1)
            jmag = INT(x)

        call mxCopyPtrToReal8(mxGetPr(PRHS(2)),jf_r,50)
        do j = 1, 50
            if (jf_r(j) .ne. 0) then
                jf(j) = .true.
            else
                jf(j) = .false.
            end if
        end do

        call mxCopyPtrToReal8(mxGetPr(prhs(3)),alati,1)
        call mxCopyPtrToReal8(mxGetPr(prhs(4)),along,1)
        call mxCopyPtrToReal8(mxGetPr(prhs(5)),x,1)
            iyyyy = INT(x)
        call mxCopyPtrToReal8(mxGetPr(prhs(6)),x,1)
            mmdd = INT(x)
        call mxCopyPtrToReal8(mxGetPr(prhs(7)),x,1)
            iut = INT(x)
        call mxCopyPtrToReal8(mxGetPr(prhs(8)),dhour,1)
        call mxCopyPtrToReal8(mxGetPr(prhs(9)),x,1) 
            height = INT(x)
        call mxCopyPtrToReal8(mxGetPr(prhs(10)),h_tec_max,1)
        call mxCopyPtrToReal8(mxGetPr(prhs(11)),x,1)
            ivar = INT(x)
        call mxCopyPtrToReal8(mxGetPr(prhs(12)),vbeg,1)
        call mxCopyPtrToReal8(mxGetPr(prhs(13)),vend,1)
        call mxCopyPtrToReal8(mxGetPr(prhs(14)),vstp,1) 
        call mxCopyPtrToReal8(mxGetPr(prhs(15)),oarr_r,100)

c =======================================================================
c Programs using subroutine IRI_SUB need to include (see IRITEST.FOR):
c
        call read_ig_rz
        call readapf107
c =======================================================================
c Programs using subroutine IRI_WEB need to include (see IRITEST.FOR):
c
c        do i=1,100
c          oarr(i,1)=-1.0
c        end do
c =======================================================================
        do j = 1, 100
           oarr(j,1) = oarr_r(j)
        end do

C =======================================================================
C CREATE MATRIXES FOR RETURN ARGUMENTS      
C =======================================================================
C OUTPUT:  OUTF(1:20,1:1000) - REAL DIMENSION
C -----------------------------------------------------------------------
        plhs(1) = mxCreateNumericMatrix(20, 1000,  
     &                  mxClassIDFromClassName('single'), 0)
C =======================================================================
C OUTPUT:  OARR(1:100,1:1000) - REAL DIMENSION
C -----------------------------------------------------------------------
        plhs(2) = mxCreateNumericMatrix(100, 1000,       
     &                  mxClassIDFromClassName('single'), 0)

C =======================================================================
C ASSIGN POINTERS TO THE VARIOUS PARAMETERS
C =======================================================================
        P_OUTF = MXGETPR(PLHS(1))
        P_OARR = MXGETPR(PLHS(2))
C =======================================================================
      call iri_web(jmag, jf, real(alati), real(along),
     &             iyyyy, mmdd, iut, real(dhour),
     &             height, real(h_tec_max), ivar,
     &             real(vbeg), real(vend), real(vstp), outf, oarr)

C
C COPY OUTPUT WHICH IS STORED IN LOCAL ARRAY TO MATRIX OUTPUT
        CALL MXCOPYREAL4TOPTR(OUTF, P_OUTF, 20*1000)
        CALL MXCOPYREAL4TOPTR(OARR, P_OARR, 100*1000)

      return
      end

C-----------------------------------------------------------------------