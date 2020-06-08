C     Calls the fortran program bsplpp that converts
C     B-form coeff to pp form.
C     simplified version with no repeated knots
C     coef = mexbsplpp(knots,fknots)
C     knots(l+7), fknots(l+3)

C     The gateway routine
      subroutine mexFunction(nlhs, plhs, nrhs, prhs)

      integer mxGetM, mxGetN, mxGetPr
      integer mxIsNumeric, mxCreateDoubleMatrix
      integer plhs(*), prhs(*)
      integer knots_pr, fknots_pr, coef_pr
      integer nlhs, nrhs
      integer m, n, size1, size2, k, l,sizecoef
      real*8 knots(1000), fknots(1000), coef(4,1000),scrtch(4,4)
      real*8 break(1000),bcoef(1000),mxGetScalar

C     Check for proper number of arguments. 
      if(nrhs .ne. 2) then
         call mexErrMsgTxt('Two vectors input required.')
      elseif(nlhs .ne. 1) then
         call mexErrMsgTxt('One output required.')
      endif

C     Get the sizes of input arrays.
      m = mxGetM(prhs(1))
      n = mxGetN(prhs(1))
      size1 = m*n
      m = mxGetM(prhs(2))
      n = mxGetN(prhs(2))
      size2 = m*n
      n = size2
      l = n - 3

C     Column * row should be smaller than 1000.
      if(size1.gt.1000 .or. size2.gt.1000) then
         call mexErrMsgTxt('Row * column must be <= 1000.')
      endif

C     Check to ensure the input is a number.
      if(mxIsNumeric(prhs(1)) .eq. 0) then
         call mexErrMsgTxt('Input must be a number.')
      endif

C     Create matrix for the return argument.
      plhs(1) = mxCreateDoubleMatrix(4, l , 0)
      coef_pr = mxGetPr(plhs(1))

C     Copy input vectors to knots and fknots
      knots_pr = mxGetPr(prhs(1))
      fknots_pr = mxGetPr(prhs(2))
      call mxCopyPtrToReal8(knots_pr, knots, size1)
      call mxCopyPtrToReal8(fknots_pr, fknots, size2)

C     Call the computational subroutine.
      k = 4
      call bsplpp(knots,fknots,size2,k,scrtch,break,coef,l)

C     Load the data into coef_pr, which is the output to MATLAB.
      sizecoef = 4*(n-3)
      call mxCopyReal8ToPtr(coef, coef_pr, sizecoef)     
      return
      end

