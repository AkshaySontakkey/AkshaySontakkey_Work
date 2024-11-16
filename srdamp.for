      subroutine vdload (
C Read only (unmodifiable)variables -
     1 nBlock, ndim, stepTime, totalTime,
     2 amplitude, curCoords, velocity, dirCos, jltyp, sname,
C Write only (modifiable) variable -
     1 value )
C
      include 'vaba_param.inc'
C
      dimension curCoords(nBlock,ndim), velocity(nBlock,ndim),
     1  dirCos(nBlock,ndim,ndim), value(nBlock)
      character*80 sname
C
      PARAMETER (PI=3.141593)
c     initiating current cords	  
      PARAMETER (R0=1.0E-3, IPEAK=40.0E3, ALPHA=0.0932, BETA=10900.D0, TM=10E-6, TR=40E-6)	
	  OPEN(unit=18, file='G:\A2\Paper3\explicit\srdata.dat', status='UNKNOWN') 
 
C     ARC CHANNEL ORIGINE
      X1=0.0D0
      Y1=0.0D0
      JLTYP=0
c	  Magnetic Permeability (H/m)
      perm=PI*4.0D-7
c     RADIUS EXPANDS DURING TR
      R=R0+ALPHA*((IPEAK)**(1.0/3))*(stepTime**(1.0/2))
C	  RADIUS IN MM
      R=R*1000
C
      IF (stepTime<=TM) THEN
		  I=IPEAK*(stepTime/TM)
	      p=(perm*I**2)/(8*(3.141593*R)**2)		  
          do km = 1, nBlock
              X=curCoords(km,1)
              Y=curCoords(km,2)
		      rsqre = X**2 + Y**2 
		      rinv = (1/rsqre)**half
              s = ceil( floor( R * rinv ) / ( floor(( R * rinv ) + 1 ))) 
              value(km)=p*s
              write (18,*) value			  
          end do		  
      ELSE
		  I=IPEAK*EXP(-BETA*(stepTime-TM))
	      p=(perm*I**2)/(8*(3.141593*R)**2)		  
          do km = 1, nBlock
              X=curCoords(km,1)
              Y=curCoords(km,2)
		      rsqre = X**2 + Y**2 
		      rinv = (1/rsqre)**half
              s = ceil( floor( R * rinv ) / ( floor(( R * rinv ) + 1 ))) 
              value(km)=p*s			  
          end do		  
      END IF
c      write (18,*) value	  
      return
      end
c
      subroutine vumat(
C Read only -
     1  nblock, ndir, nshr, nstatev, nfieldv, nprops, jInfoArray,
     2  stepTime, totalTime, dtArray, cmname, coordMp, charLength,
     3  props, density, strainInc, relSpinInc,
     4  tempOld, stretchOld, defgradOld, fieldOld,
     3  stressOld, stateOld, enerInternOld, enerInelasOld,
     6  tempNew, stretchNew, defgradNew, fieldNew,
C Write only -
     5  stressNew, stateNew, enerInternNew, enerInelasNew )
C
      include 'vaba_param.inc'
c 
c 3D Orthotropic Elasticity with Hashin 3d Failure criterion 
c 
c The state variables are stored as: 
c    state(*,1)   = material point status 
c    state(*,2:7) = damping stresses 
c    state(*,3)   = temperature 
c 
c User defined material properties are stored as 
c  * First line: 
c     props(1) --> Young's modulus in 1-direction, E1 
c     props(2) --> Young's modulus in 2-direction, E2,E3 
c     props(3) --> Poisson's ratio, nu12, nu13 
c     props(4) --> Poisson's ratio, nu23 
c     props(5) --> Shear modulus, G12, G13 
c     props(6) --> Shear modulus, G23 
c     props(7) --> alpha11
c     props(8) --> alpha22, alpha33 
c 
c  * Second line: 
c     props(9)  --> "not used" 
c     props(10) --> beta damping parameter 
c     props(11) --> Young's modulus in 3-direction, E1 at T(ablation) 
c     props(12) --> Young's modulus in 2-direction, E2, E3 at T(200 + off) 
c     props(13) --> Young's modulus in 2-direction, E2, E3 at T(260 + off) 
c     props(14) --> Young's modulus in 2-direction, E2, E3 at T(600 + off) 
c     props(15) --> Young's modulus in 2-direction, E2, E3 at T(ablation) 
c     props(16) --> "not used" 
c 
c  * Third line: 
c     props(17) --> Ultimate tens stress in 1-direction, sigu1t 
c     props(18) --> Ultimate comp stress in 1-direction, sigu1c 
c     props(19) --> Ultimate tens stress in 2&3-direction, sigu2t, sigu3t 
c     props(20) --> Ultimate comp stress in 2&3-direction, sigu2c, sigu3c 
c     props(21) --> Ultimate shear stress, sigu12, sigu13, sigu23  
c     props(22) --> "not used" 
c     props(23) --> Shear modulus, G12, G13 at T(200 + 0ff)
c     props(24) --> Shear modulus, G12, G13 at T(260 + 0ff) 
c 
c  * Fourth line: 
c     props(25) --> Shear modulus, G12, G13 at T(600 + 0ff)  
c     props(26) --> Shear modulus, G12, G13 at T(ablation)  
c     props(27) --> Shear modulus, G23 at T(200 + 0ff)  
c     props(28) --> Shear modulus, G23 at T(260 + 0ff) 
c     props(29) --> Shear modulus, G23 at T(600 + 0ff) 
c     props(30) --> Shear modulus, G23 at T(ablation) 
c     props(31) --> "not used" 
c     props(32) --> "not used"
c     props(33) --> "Fracture Energy GFc Fibre compression 1400 J/m² Kumar Carbone/Epoxy"
c     props(34) --> "Fracture Energy GFt Fibre tension 1400J/m² Kumar for Carbone/Epoxy"
c     props(35) --> "Fracture Energy GMc Matrix compression 200 J/m² Kumar Carbone/Epoxy"      
c     props(36) --> "Fracture Energy GMt Matrix tension 7000 J/m² Kumar Carbone/Epoxy"
c
c Declaration of Critic Fracture Energy GFt,GFc, GMt and GMc
c
      dimension props(nprops), density(nblock), coordMp(nblock,*),
     1  charLength(nblock), strainInc(nblock,ndir+nshr),
     2  relSpinInc(nblock,nshr), tempOld(nblock), 
     3  stretchOld(nblock,ndir+nshr),
     4  defgradOld(nblock,ndir+nshr+nshr),
     5  fieldOld(nblock,nfieldv), stressOld(nblock,ndir+nshr),
     6  stateOld(nblock,nstatev), enerInternOld(nblock),
     7  enerInelasOld(nblock), tempNew(nblock),
     8  stretchNew(nblock,ndir+nshr),
     8  defgradNew(nblock,ndir+nshr+nshr),
     9  fieldNew(nblock,nfieldv),
     1  stressNew(nblock,ndir+nshr), stateNew(nblock,nstatev),
     2  enerInternNew(nblock), enerInelasNew(nblock), jInfoArray(*)
C All arrays dimensioned by (*) are not used in this algorithm
C
      character*80 cmname
c  
* Temporary arrays 
      dimension eigen(maxblk*3), statusMp(nblock), dtstrain(nblock,6),
     *     alpha11(nblock), alpha22(nblock), alpha33(nblock), SR(ndir + nshr),
     *     C11(nblock), C22(nblock), C33(nblock), 
     *     C12(nblock), C13(nblock), C23(nblock),	 
     *     G12(nblock), G13(nblock), G23(nblock), curDensity(nBlock)	 
	  DOUBLE PRECISION SR2, sf1, sf2, sf3, sf4, sf5, sf6
	  real*8 SRFactor, strainIncMag, SRIncMag
*
* 
      parameter(
     *	i_svd_DmgFiberT  = 1, i_svd_DmgFiberC = 2,
     *	i_svd_DmgMatrixT = 3, i_svd_DmgMatrixC = 4,
     *	i_svd_statusMp   = 5, i_svd_Strain = 6,
     *    i_svd_enomMax = 12, i_svd_enomMin = 13,      
     *    i_svd_dampStress = 14,
     *	n_svd_required = 19 )
*
      parameter( 
     *     i_s33_Xx = 1, 
     *     i_s33_Yy = 2, 
     *     i_s33_Zz = 3, 
     *     i_s33_Xy = 4, 
     *     i_s33_Yz = 5, 
     *     i_s33_Zx = 6,
     *     i_pro_fenegfc = 33,	 
     *     i_pro_fenegft = 34,
     *     i_pro_fenegmc = 35,
     *     i_pro_fenegmt = 36) 
*
      parameter( zero = 0.d0, one = 1.d0, two = 2.d0, half = 0.5d0, pontone = 0.1d0,
     *     reftemp = 25.d0, twoh = 200.d0, twosxty = 260.d0, oneh = 100.d0,
     *     sixh = 600.d0, ablation = 3316.d0, off = 281.d0, sfmodule = 1.29d0, sfstress = 1.38d0, sfftc = 1.6d0,
     *     three = 3.d0, onpsf = 1.75d0, sffttwoc = 1.25d0, sfftt = 1.23, sffttwot = 1.2)    
* Read material properties
      xnu12 = props(3) 
      xnu13 = props(3) 
      xnu23 = props(4)      	  	 
*
c      DO i=1, (n_s33_Car)
      SR2=0.0d0
c      end do
c Compute the magnitude of strain increment
c     strainIncMag = 0.0d0
	  SRIncMag = 0.0d0
c     do k = 1, nblock
c         do j=1, 6
C             strainIncMag = strainIncMag + strainInc(k,j)**2
C             if (j .eq. 2) then 
             SRIncMag = SRIncMag + strainInc(k,2)**2
C             end if
C         end do
c      end do
	  SRIncMag = sqrt(SRIncMag)
c Compute the strain rate
c     DO i = 1, 6
      SR(2) = SRIncMag/dtArray
	  SR2 = SR(2)
c     END DO
    

* compute strain rate dependent properties in e2
c		sf1 - scale factor for E2, E3 sfmodule = 1.29
c       sf2 - scale factor for xt, xc, yt, yc sfstress = 1.38
c	    sf3 - scale factor for fracture toughness sfftc = 1.6 Quasi static = 10 high rate = 16
c		sf4 - scale factor for fracture toughness sfftt = 1.23 Quasi static = 133 high rate = 164
c	    sf5 - scale factor for fracture toughness sffttwoc = 1.25 Quasi static = 1.6 high rate = 2
c		sf6 - scale factor for fracture toughness sffttwot = 1.2 Quasi static = 0.5 high rate = 0.6
		  if (SR2 <= pontone) then 
		      sf1 = one
			  sf2 = one
			  sf3 = one
			  sf4 = one
			  sf5 = one
			  sf6 = one
		  else if (SR2 < oneh) then 
		      sf1 = ((SR2 - pontone)*(sfmodule - one)/(oneh - pontone)) + one
			  sf2 = ((SR2 - pontone)*(sfstress - one)/(oneh - pontone)) + one
			  sf3 = ((SR2 - pontone)*(sfftc - one)/(oneh - pontone)) + one
			  sf4 = ((SR2 - pontone)*(sfftt - one)/(oneh - pontone)) + one
			  sf5 = ((SR2 - pontone)*(sffttwoc - one)/(oneh - pontone)) + one
			  sf6 = ((SR2 - pontone)*(sffttwot - one)/(oneh - pontone)) + one
		  else 
		      sf1 = sfmodule
			  sf2 = sfstress
			  sf3 = sfftc
			  sf4 = sfftt
			  sf5 = sffttwoc
			  sf6 = sffttwot
		  end if
c assign properties to thermal expansion coefficient
      do k = 1, nblock
	      if (tempOld(k) .le. (twoh + off)) then
              alpha11(k) = props(7) 
              alpha22(k) = props(8)
              alpha33(k) = props(8)
			  curDensity(k) = props(31)
          else if (tempOld(k) .gt. (twoh + off) .and. tempOld(k) .le. ablation) then
              fac1=((tempOld(k)-(twoh + off))/(sixh - twoh))
			  If (fac1 .lt. zero) then fac1 = zero end if
			  If (fac1 .gt. one) then fac1 = one end if
			  fac0=one - fac1		  
              alpha11(k) = props(7)*fac0 + fac1 * props(7) * three
              alpha22(k) = props(8)*fac0 + fac1 * props(8) * onpsf
              alpha33(k) = props(8)*fac0 + fac1 * props(8) * onpsf
			  curDensity(k) = props(31)*fac0 + fac1 * props(32)
          else if (tempOld(k) .gt.  ablation) then
		      stateNew(k,i_svd_statusMp) = zero
              alpha11(k) = props(7) * three
              alpha22(k) = props(8) * onpsfs
              alpha33(k) = props(8) * onpsf
			  curDensity(k) = 1.11d-12			  
          end if
*		  
	  end do
*
c Assign temperature dependent material and heat rate, strain rate scale factor 
      do k = 1, nblock  
          if (tempOld(k) .le. (twoh + off)) then
              E1 = props(1)			  
              E2 = props(2) * sf1
              E3 = props(2) * sf1
              G12(k) = props(5)
              G13(k) = props(5)
              G23(k) = props(6)
          else if (tempOld(k) .gt. (twoh + off) .and. tempOld(k) 
     *            .le. (twosxty + off)) then
              fac1=((tempOld(k)-(twoh + off))/(twosxty - twoh))
			  If (fac1 .lt. zero) then fac1 = zero end if
			  If (fac1 .gt. one) then fac1 = one end if
			  fac0 = one - fac1
              E1 = props(1)		 
              E2 = (props(12)*fac0 + fac1 * props(13)) * sf1
              E3 = (props(12)*fac0 + fac1 * props(13)) * sf1
              G12(k) = props(23)*fac0 + fac1 * props(24) 
              G13(k) = props(23)*fac0 + fac1 * props(24)
              G23(k) = props(27)*fac0 + fac1 * props(28)	
          else if (tempOld(k) .gt. (twosxty + off) .and. tempOld(k) 
     *            .le. (sixh + off)) then
              fac1=((tempOld(k)-(twosxty + off))/(sixh - twosxty))
			  If (fac1 .lt. zero) then fac1 = zero end if
			  If (fac1 .gt. one) then fac1 = one end if
			  fac0=one - fac1
              E1 = props(1)		 
              E2 = (props(13)*fac0 + fac1 * props(14)) * sf1
              E3 = (props(13)*fac0 + fac1 * props(14)) * sf1
              G12(k) = props(24)*fac0 + fac1 * props(25) 
              G13(k) = props(24)*fac0 + fac1 * props(25) 
              G23(k) = props(28)*fac0 + fac1 * props(29)             
          else if (tempOld(k) .gt. (sixh + off) .and. tempOld(k) .le. ablation) then
              E1 = props(1)		 
              E2 = props(14) * sf1
              E3 = props(14) * sf1
              G12(k) = props(25) 
              G13(k) = props(25) 
              G23(k) = props(29)
          else if (tempOld(k) .gt. ablation) then
              E1 = props(11)		 
              E2 = props(15) * sf1
              E3 = props(15) * sf1
              G12(k) = props(26) 
              G13(k) = props(26) 
              G23(k) = props(30)			  
          end if
* 
          xnu21 = xnu12 * E2 / E1 
          xnu31 = xnu13 * E3 / E1 
          xnu32 = xnu23 * E3 / E2 
* 		 
* Compute terms of stiffness matrix 
          gg = one / ( one - xnu12*xnu21 - xnu23*xnu32 - xnu31*xnu13 
     *     - two*xnu21*xnu32*xnu13 ) 
          C11(k)  = E1 * ( one - xnu23*xnu32 ) * gg 
          C22(k)  = E2 * ( one - xnu13*xnu31 ) * gg 
          C33(k)  = E3 * ( one - xnu12*xnu21 ) * gg 
          C12(k)  = E1 * ( xnu21 + xnu31*xnu23 ) * gg
          C13(k)  = E1 * ( xnu31 + xnu21*xnu32 ) * gg 
          C23(k)  = E2 * ( xnu32 + xnu12*xnu31 ) * gg 
      end do
		  
* 
      f1t = props(17) * sf2
      f1c = props(18) * sf2
      f2t = props(19) * sf2
      f2c = props(20) * sf2
      f3t = props(19) * sf2
      f3c = props(20) * sf2
      f12 = props(21) * sf2
      f13 = props(21) * sf2
      f23 = props(21) * sf2
c Declaration of Critic Fracture Energy GFt,GFc, GMt and GMc
      fenegfc = props(i_pro_fenegfc) * sf3
      fenegft = props(i_pro_fenegft) * sf4
      fenegmc = props(i_pro_fenegmc) * sf5
      fenegmt = props(i_pro_fenegmt) * sf6
c      print*,"charlength_before_calls=",charlength	  
* 
      beta = props(10)
c      print*,"stretchNew1111",stretchNew	  
* 
* Assume purely elastic material at the beginning of the analysis
*       
      if ( totalTime .eq. zero ) then	  
         if (nstatev .lt. n_svd_Required) then 
            call xplb_abqerr(-2,'Subroutine VUMAT requires the '// 
     *           'specification of %I state variables. Check the '// 
     *           'definition of *DEPVAR in the input file.', 
     *           n_svd_Required,zero,' ') 
            call xplb_exit 
         end if 
         sf1 = one
		 sf2 = one
		 sf3 = one
		 sf4 = one
		 sf5 = one
		 sf6 = one         
		 call OrthoEla3dExp ( nblock, 
     *        stateOld(1,i_svd_DmgFiberT), 
     *        stateOld(1,i_svd_DmgFiberC), 
     *        stateOld(1,i_svd_DmgMatrixT), 
     *        stateOld(1,i_svd_DmgMatrixC), 
     *        C11, C22, C33, C12, C23, C13, G12, G23, G13, 
     *        strainInc, 
     *        stressNew )
	     call thermalstrain( nblock, tempOld, alpha11, alpha22, 
     *        alpha33, stateOld(1,i_svd_strain) )	 
         return		 
      end if 	  
* 
      do k=1, nBlock
          dtstrain(k,i_s33_Xx)=	 alpha11(k)*(tempNew(k)-tempOld(k))
          dtstrain(k,i_s33_Yy)=	 alpha22(k)*(tempNew(k)-tempOld(k))
          dtstrain(k,i_s33_Zz)=	 alpha33(k)*(tempNew(k)-tempOld(k))
      end do		  
*  Update total elastic strain 
      call strainUpdate ( nblock, strainInc, 
     *     stateOld(1,i_svd_strain), stateNew(1,i_svd_strain), dtstrain ) 
* 
* Stress update 
      call OrthoEla3dExp ( nblock, 
     *     stateOld(1,i_svd_DmgFiberT), 
     *     stateOld(1,i_svd_DmgFiberC), 
     *     stateOld(1,i_svd_DmgMatrixT), 
     *     stateOld(1,i_svd_DmgMatrixC), 
     *     C11, C22, C33, C12, C23, C13, G12, G23, G13, 
     *     stateNew(1,i_svd_strain), 
     *     stressNew ) 
* 
* Failure evaluation 
* 
      call copyr ( nblock, 
     *     stateOld(1,i_svd_DmgFiberT), stateNew(1,i_svd_DmgFiberT) ) 
      call copyr ( nblock, 
     *     stateOld(1,i_svd_DmgFiberC), stateNew(1,i_svd_DmgFiberC) ) 
      call copyr ( nblock, 
     *     stateOld(1,i_svd_DmgMatrixT), stateNew(1,i_svd_DmgMatrixT) ) 
      call copyr ( nblock, 
     *     stateOld(1,i_svd_DmgMatrixC), stateNew(1,i_svd_DmgMatrixC) ) 
      call copyr ( nblock,
     1	stateOld(1,i_svd_enomMax), stateNew(1,i_svd_enomMax) )
      call copyr ( nblock,
     1	stateOld(1,i_svd_enomMin), stateNew(1,i_svd_enomMin) )
      nDmg = 0 
      call eig33Anal ( nblock, stretchNew, eigen ) 
      call Hashin3d ( nblock, nDmg,
     *     f1t, f2t, f3t, f1c, f2c, f3c, f12, f23, f13,
     *     stateNew(1,i_svd_DmgFiberT), stateNew(1,i_svd_DmgFiberC),
     *     stateNew(1,i_svd_DmgMatrixT), stateNew(1,i_svd_DmgMatrixC),
     *     stateNew(1,i_svd_statusMp),
     *	   stressNew,eigen, stateNew(1,i_svd_strain),
     *     charLength, fenegft, fenegfc, fenegmt, fenegmc,
     *     E1, E2, E3,stateNew(1,i_svd_enomMax),
     *     stateNew(1,i_svd_enomMin) )
*     -- Recompute stresses if new Damage is occurring 
      if ( nDmg .gt. 0 ) then 
         call OrthoEla3dExp ( nblock, 
     *        stateNew(1,i_svd_DmgFiberT), 
     *        stateNew(1,i_svd_DmgFiberC), 
     *        stateNew(1,i_svd_DmgMatrixT), 
     *        stateNew(1,i_svd_DmgMatrixC), 
     *        C11, C22, C33, C12, C23, C13, G12, G23, G13, 
     *        stateNew(1,i_svd_strain), 
     *        stressNew ) 
      end if 
* 
* Beta damping 
      if ( beta .gt. zero ) then 
         call betaDamping3d ( nblock, 
     *        beta, dtArray, strainInc, 
     *        stressOld, stressNew, 
     *        stateNew(1,i_svd_statusMp), 
     *        stateOld(1,i_svd_dampStress), 
     *        stateNew(1,i_svd_dampStress) ) 
      end if 
* 
* Integrate the internal specific energy (per unit mass) 
* 
      call EnergyInternal3d ( nblock, stressOld, stressNew,
     *   strainInc, curDensity, enerInternOld, enerInternNew) 
* 
      return 
      end 
c
************************************************************ 
*   strainUpdate_initial: Update thermal strain                      * 
************************************************************ 
      subroutine thermalstrain ( nblock, tempOld, alpha11, alpha22, 
     *     alpha33, strain ) 
* 
      include 'vaba_param.inc' 
* 
      parameter( reftemp = 25.d0, mone= -1.d0 )
      parameter(
     *     i_s33_Xx = 1, 
     *     i_s33_Yy = 2, 
     *     i_s33_Zz = 3, 
     *     i_s33_Xy = 4, 
     *     i_s33_Yz = 5, 
     *     i_s33_Zx = 6, 
     *     n_s33_Car = 6 ) 
* 
      dimension strain(nblock,n_s33_Car),
     *     tempOld(nblock), alpha11(nblock), alpha22(nblock), alpha33(nblock)	 
*
      do k=1, nBlock 
          strain(k,i_s33_Xx)= mone*alpha11(k)*(tempOld(k)-reftemp) 
          strain(k,i_s33_Yy)= mone*alpha22(k)*(tempOld(k)-reftemp) 
          strain(k,i_s33_Zz)= mone*alpha33(k)*(tempOld(k)-reftemp)		  
      end do 
* 
      return 
      end 
************************************************************ 
*   OrthoEla3dExp: Orthotropic elasticity - 3d             * 
************************************************************ 
      subroutine OrthoEla3dExp ( nblock, 
     *     dmgFiberT, dmgFiberC, dmgMatrixT, dmgMatrixC, 
     *     C11, C22, C33, C12, C23, C13, G12, G23, G13, 
     *     strain, stress ) 
* 
      include 'vaba_param.inc' 

*  Orthotropic elasticity, 3D case - 
* 
      parameter( zero = 0.d0, one = 1.d0, two = 2.d0) 
      parameter( 
     *     i_s33_Xx = 1, 
     *     i_s33_Yy = 2, 
     *     i_s33_Zz = 3, 
     *     i_s33_Xy = 4, 
     *     i_s33_Yz = 5, 
     *     i_s33_Zx = 6, 
     *     n_s33_Car = 6 ) 
* 
      dimension  strain(nblock,n_s33_Car), 
     *     dmgFiberT(nblock), dmgFiberC(nblock), 
     *     dmgMatrixT(nblock), dmgMatrixC(nblock), 
     *     stress(nblock,n_s33_Car) 
      dimension C11(nblock), C22(nblock), C33(nblock), 
     *     C12(nblock), C13(nblock), C23(nblock),	 
     *     G12(nblock), G13(nblock), G23(nblock)	 
*     -- shear fraction in matrix tension and compression mode 
      parameter( smt = 0.9d0, smc = 0.5d0 ) 
* 
      do k = 1, nblock 
*     -- Compute damaged stiffness 
          dft = dmgFiberT(k) 
          dfc = dmgFiberC(k) 
          dmt = dmgMatrixT(k) 
          dmc = dmgMatrixC(k) 
          df = one - ( one - dft ) * ( one - dfc ) 
* 
          dC11 = ( one - df ) * C11(k) 
          dC22 = ( one - df ) * ( one - dmt ) * ( one - dmc ) * C22(k) 
          dC33 = ( one - df ) * ( one - dmt ) * ( one - dmc ) * C33(k) 
          dC12 = ( one - df ) * ( one - dmt ) * ( one - dmc ) * C12(k) 
          dC23 = ( one - df ) * ( one - dmt ) * ( one - dmc ) * C23(k) 
          dC13 = ( one - df ) * ( one - dmt ) * ( one - dmc ) * C13(k) 
          dG12 = ( one - df ) 
     *        * ( one - smt*dmt ) * ( one - smc*dmc ) * G12(k) 
          dG23 = ( one - df ) 
     *        * ( one - smt*dmt ) * ( one - smc*dmc ) * G23(k) 
          dG13 = ( one - df ) 
     *        * ( one - smt*dmt ) * ( one - smc*dmc ) * G13(k) 
*     -- Stress update 
          stress(k,i_s33_Xx) = dC11 * strain(k,i_s33_Xx) 
     *        + dC12 * strain(k,i_s33_Yy) 
     *        + dC13 * strain(k,i_s33_Zz) 
          stress(k,i_s33_Yy) = dC12 * strain(k,i_s33_Xx) 
     *        + dC22 * strain(k,i_s33_Yy) 
     *        + dC23 * strain(k,i_s33_Zz) 
          stress(k,i_s33_Zz) = dC13 * strain(k,i_s33_Xx) 
     *        + dC23 * strain(k,i_s33_Yy) 
     *        + dC33 * strain(k,i_s33_Zz) 
          stress(k,i_s33_Xy) = two * dG12 * strain(k,i_s33_Xy) 
          stress(k,i_s33_Yz) = two * dG23 * strain(k,i_s33_Yz) 
          stress(k,i_s33_Zx) = two * dG13 * strain(k,i_s33_Zx) 
      end do 
*     
      return 
      end 
************************************************************ 
*   strainUpdate: Update total strain                      * 
************************************************************ 
      subroutine strainUpdate ( nblock, 
     *     strainInc, strainOld, strainNew, dtstrain ) 
* 
      include 'vaba_param.inc' 
* 
      parameter( 
     *     i_s33_Xx = 1, 
     *     i_s33_Yy = 2, 
     *     i_s33_Zz = 3, 
     *     i_s33_Xy = 4, 
     *     i_s33_Yz = 5, 
     *     i_s33_Zx = 6, 
     *     n_s33_Car = 6 ) 
* 
      dimension strainInc(nblock,n_s33_Car), dtstrain(nblock,n_s33_Car), 
     *     strainOld(nblock,n_s33_Car), 
     *     strainNew(nblock,n_s33_Car) 
* 
      do k = 1, nblock 
          strainNew(k,i_s33_Xx)= strainOld(k,i_s33_Xx) 
     *                        + strainInc(k,i_s33_Xx) + dtstrain(nblock,i_s33_Xx)
          strainNew(k,i_s33_Yy)= strainOld(k,i_s33_Yy) 
     *                        + strainInc(k,i_s33_Yy) + dtstrain(nblock,i_s33_Yy)
          strainNew(k,i_s33_Zz)= strainOld(k,i_s33_Zz) 
     *                        + strainInc(k,i_s33_Zz) + dtstrain(nblock,i_s33_Zz)
          strainNew(k,i_s33_Xy)= strainOld(k,i_s33_Xy) 
     *                        + strainInc(k,i_s33_Xy) 
          strainNew(k,i_s33_Yz)= strainOld(k,i_s33_Yz) 
     *                        + strainInc(k,i_s33_Yz) 
          strainNew(k,i_s33_Zx)= strainOld(k,i_s33_Zx) 
     *                        + strainInc(k,i_s33_Zx) 
      end do 
* 
      return 
      end 
************************************************************ 
*   Hashin3d with Modified Puck: Evaluate Hashin/ Puck 3d failure  * 
*   criterion for fiber, Puck for matrix                   * 
************************************************************ 
c      
      subroutine Hashin3d ( nblock, nDmg,
     *    f1t, f2t, f3t, f1c, f2c, f3c, f12, f23, f13,
     *    dmgFiberT, dmgFiberC, dmgMatrixT, dmgMatrixC,
     *    statusMp, stress,eigen, strain,charLength,
     *    fenegft, fenegfc, fenegmt, fenegmc,
     *    E1, E2, E3, enomMax, enomMin ) 
c -------------
      include 'vaba_param.inc'
c -------------
      parameter( zero  = 0.d0, one = 1.d0, half = 0.5d0,
     *           three = 3.d0, two = 2.d0, X = 0.d0 )
c -------------
      parameter(
     1 i_s33_Xx = 1,
     2 i_s33_Yy = 2,
     3 i_s33_Zz = 3,
     4 i_s33_Xy = 4,
     5 i_s33_Yz = 5,
     6 i_s33_Zx = 6,
     7 n_s33_Car = 6 )
      parameter(i_v3d_X=1,i_v3d_Y=2,i_v3d_Z=3 )
      parameter(n_v3d_Car=3 )
c
c* NB: The default limiting strain values are eMax 1.00d0 maximum and eMin -0.8d0 minimum, and are*
*	defined within the VUMAT; these limits are in place to minimise
computational*******************
*	complications from heavily distorted elements rather than to describe the physical**************
*	behaviour of the laminate***********************************************************************
      parameter ( eMax = 1.00d0, eMin = -0.8d0 )
c -------------
      dimension dmgFiberT(nblock), dmgFiberC(nblock),
     *	      dmgMatrixT(nblock), dmgMatrixC(nblock),
c --------      
     *	statusMp(nblock), strain(nblock,n_s33_Car),
     *	stress(nblock,n_s33_Car), eigen(nblock,n_s33_Car),DmgMax(nblock),
     *    charLength(nblock),
c-----
     *     enomMax(nblock), enomMin(nblock)
c---------
c      print*,"charleng_Daroon_Sub_Hashin=",charlength,nblock
c      print*,"dis0ft111",dis0ft
      f1tInv = zero
      f2tInv = zero
      f3tInv = zero
      f1cInv = zero
      f2cInv = zero
      f3cInv = zero
      f12Inv = zero
      f23Inv = zero
      f13Inv = zero
c
      if ( f1t .gt. zero ) f1tInv = one / f1t
      if ( f2t .gt. zero ) f2tInv = one / f2t
      if ( f3t .gt. zero ) f3tInv = one / f3t
      if ( f1c .gt. zero ) f1cInv = one / f1c
      if ( f2c .gt. zero ) f2cInv = one / f2c
      if ( f3c .gt. zero ) f3cInv = one / f3c
      if ( f12 .gt. zero ) f12Inv = one / f12
      if ( f23 .gt. zero ) f23Inv = one / f23
      if ( f13 .gt. zero ) f13Inv = one / f13
c
      do k = 1, nblock
      if ( statusMp(k) .eq. one ) then
c          print*,"statusmp =", statusMp(k)
c
      ldmg = 0
c
      s11 = stress(k,i_s33_Xx)
      s22 = stress(k,i_s33_Yy)
      s33 = stress(k,i_s33_Zz)
      s12 = stress(k,i_s33_Xy)
      s23 = stress(k,i_s33_Yz)
      s13 = stress(k,i_s33_Zx)
      STRAIN11 = STRAIN(K,I_S33_XX)
      STRAIN22 = STRAIN(K,I_S33_YY)
      STRAIN33 = STRAIN(K,I_S33_ZZ)
      STRAIN12 = STRAIN(K,I_S33_XY)
      STRAIN23 = STRAIN(K,I_S33_YZ)
      STRAIN13 = STRAIN(K,I_S33_ZX)
c      charL    = charLength(k)
c      print*, "Charlength3=",charl,charlength(k)
c============================================
c               Fiber Mode
c============================================
c Evaluate Fiber modes
      if ( s11 .gt. zero ) then
c -- Tensile Fiber Mode
      rft = (s11*f1tInv)**2 + (s12*f12Inv)**2 + (s13*f13Inv)**2
      if ( rft .ge. one ) then
c
**************************Evaluate Fibre Modes******************************************************
*Evaluate Fibre damage Tension: Mode 1*
*******Compute Equivalent Displacement in the Second Part of the Curve Stress-Displacement (diseqft)*
       diseqft = charLength(k) * sqrt((STRAIN11**2 +
     1 STRAIN12**2 + STRAIN13**2))
*******Compute Equivalent Stress in the second part of the curve Force-Displacement (streqft)*
       streqft = (CharLength(k) * (s11*
     1          STRAIN11+ s12 * STRAIN12
     2          + s13 * STRAIN13))/diseqft
** compute initial stress at maximum load following Xin 2015(stressinatft) sigma(init)=sigma(i)* f(i)
       stressinatft = streqft / sqrt(rft)
** compute initial displacement at maximum load following Xin 2015(disinatft) delta(init)=delta(i)/f(i)
       disinatft = diseqft / sqrt(rft)
***********Compute Failure Displacement Based on Fracture Energy for Mode 1 (disfailft-disinatft=2(G1)/stressinatft)
       disfailft= (2* fenegft/stressinatft)
***********Compute Damage Variable for Fibre Tension (mode1 = DmgFiberT)
       DmgFiberT(k) = (disfailft*(diseqft-disinatft))/
     1    (diseqft*(disfailft-disinatft))
      endif
 
      if (DmgFiberT(K) .ge. one) then
      lDmg = 1
      DmgFiberT(k) = one
 
      endif
      
      elseif ( s11 .lt. zero ) then
 
c**************************************************Compressive Fibre mode2: (abs(S11)+<S33>)**2/(Xc)**2
      rfc =(abs(s11)+ZZMACCAULEY(s33))**2 * (f1cInv**2)
      if ( rfc .ge. one ) then
**************************************************************************************
********************************************Compute Equivalent Displacement in the Second Part of the Curve Force-Displacement (diseqfc)**********
      diseqfc = CharLength(k) * sqrt ((-1*STRAIN11)**2)
**************Compute Equivalent Stress in the Second Part of the Curve Force-Displacement (streqfc)****************
      streqfc = CharLength(k) * ((-1*(STRAIN11))*s11)/(diseqfc)
************compute initial stress at maximum load following Xin 2015 (stressinatfc)sigma(init)=sigma(i)* f(i)*****
      stressinatfc = streqfc / sqrt(rfc)
************compute initial displacement at maximum load following Xin 2015 (disinatfc) delta(init)=delta(i)/f(i)***
      disinatfc = diseqfc / sqrt(rfc)
***********compute Failure Displacement Based on Fracture Energy for Mode 2 (disfailfc-disinatfc=2(G2)/stressinatfc)
      disfailfc = (-2 * fenegfc/stressinatfc)     !(-)ezafe Shode
************compute Damage variable for mode 2 fibre compression (DmgFiberC)
      DmgFiberC(k) = (disfailfc*(diseqfc-disinatfc))/
     1 (diseqfc*(disfailfc-disinatfc))
c
**************************************************************************************
******************************
      endif
      if (DmgFiberC(K) .ge. one) then
      lDmg = 1
      DmgFiberC(k) = one
      endif
      endif 

c Evaluate Matrix Modes
      if ( ( s22 + s33 ) .gt. zero ) then
c Tensile Matrix mode according to Hashin 3D matrix tension
c      rmt = ( s22 * f2tInv )**2
c     1	+ ( s12 * f12Inv )**2
c     2	+ ( s23 * f23Inv )**2
cHashin           rmt = ((s22 +s33) * f2tInv )**2
c     1	+ ( s23**2 - s22 * s33) * f23Inv + ( s12 * f12Inv )**2
c     2	+ ( s13 * f13Inv )**2
c* unifiberCodeFromeAbaqusDocumentation    -- Tensile Matrix mode
                  rmt = ( s11 * half * f1tInv )**2 
     *                + ( ( s22 + s33 )**2 * abs(f2tInv * f2cInv) )
     *                + ( s12 * f12Inv )**2	 
     *                + ( ( s22 + s33 ) * (f2tInv + f2cInv) )
      if ( rmt .ge. one ) then
**************************************************************************************
******************
******Compute Displacement Equivalent in the Second Part of the Curve Force- Displacement (diseqmt)
      diseqmt = CharLength(k) * sqrt(STRAIN22**2 +
     1 STRAIN12**2+ STRAIN23**2)
******Compute Equivalent Stress in the Second Part of the Curve Force-Displacement (streqmt)
      streqmt = (CharLength(k) *(s22 *
     1 STRAIN22+ s12 * STRAIN12
     2 + s23 * STRAIN23))/diseqmt
** compute initial stress at maximum load following Xin 2015 (stressinatmt)
      stressinatmt =streqmt/sqrt(rmt)
** compute initial displacement at maximum load following Xin 2015 (disinatmt)
      disinatmt = diseqmt/sqrt(rmt)
** compute failure displacement based on fracture energy for mode 3 (disfailmt- disinatmt=2(G3)/stressinatmt
      disfailmt= (2* fenegmt/stressinatmt)
** compute Damage variable for mode 2 fiber compression (DmgMatrixT)
      DmgMatrixT(k) = (disfailmt*(diseqmt - disinatmt))/         
     1 (diseqmt*(disfailmt-disinatmt))
**************************************************************************************
*******************
      endif
      if (DmgMatrixT(K) .ge. one) then
      lDmg = 1
       DmgMatrixT(k) = one
      endif
c
      elseif ( ( s22 + s33 ) .lt. zero ) then
c  Compressive Matrix Mode Hashin 3D
c      rmc =(( abs(s22)**2)*((f2cInv) **2))
c     1	+ ( s12 * f12Inv )**2
c     2	+ ( s23 * f23Inv )**2
c---*     -- Compressive Matrix Mode
c Hashin     rmc = ((s22 + s33) * f23Inv/2)**2 + (s22 + s33) * f2cInv 
c     1    * (( f23Inv / (2 * f2cInv))**2 - 1 )
c     2    + ( s23**2 - s22 * s33) * f23Inv**2 
c     3	+ ( s12 * f12Inv )**2
c     4	+ ( s23 * f23Inv )**2
c unifiberCodeFromeAbaqusDocumentation     -- Compressive Matrix Mode
                  rmc = ( s11 * half * f1tInv )**2 
     *                + ( ( s22 + s33 )**2 * abs(f2tInv * f2cInv) ) 
     *                + ( s12 * f12Inv )**2 
     *                + ( ( s22 + s33 ) * (f2tInv + f2cInv) )
      if ( rmc .ge. one ) then
**************************************************************************************
**********************
** compute displacement equivalent in the second part of the curve Force-Displacement (diseqmc)
      diseqmc = CharLength(k) * sqrt((-1*STRAIN22)**2
     1	    + (STRAIN12)**2 + (STRAIN23)**2)
** compute stress equivalent in the second part of the curve Force-Displacement (streqmc)
      streqmc = CharLength(k)*(((-1*STRAIN22) * s22)
     1	    +  (s12 * STRAIN12) 
     2	    +  (s23 * STRAIN23))/diseqmc
** compute initial stress at maximum load following Xin 2015 (stressinatmc)
      stressinatmc = streqmc / sqrt(rmc)
** compute initial displacement at maximum load following Xin 2015 (disinatmc)
      disinatmc = diseqmc /sqrt(rmc)
** compute failure displacement based on fracture energy for mode 4: disfailmc- disinatmc=2(G4)/stressinatmc
      disfailmc = (2* fenegmc/stressinatmc)
** compute Damage variable for mode 4 Matrix compression (DmgMatrixC)
      DmgMatrixC(k) = (disfailmc * (diseqmc-disinatmc))
     1           /  (diseqmc * (disfailmc-disinatmc))
**************************************************************************************
***********************
      endif
c
      end if
c
      if (DmgMatrixC(K) .ge. one) then
      lDmg = 1
      DmgMatrixC(k) = one
      endif
      eigMax  = max(eigen(k,i_v3d_X),eigen(k,i_v3d_Y),eigen(k,i_v3d_Z))
      eigMin  = min(eigen(k,i_v3d_X),eigen(k,i_v3d_Y),eigen(k,i_v3d_Z))
      enomMax(k) = eigMax - one
      enomMin(k) = eigMin - one
c Condition	ajouter par Mahrez Septembre 2015
c		DmgMax(k)= max(DmgFiberC(k),DmgFiberT(k),
c	1		DmgMatrixT(k),DmgMatrixC(k))
c	    if (DmgMax(k) .ge. one) then
c		    lDmg = 1
c	    Dmg =Max= one
c	    end if
 
          if ( enomMax(k) .gt. eMax .or. 
     1        enomMin(k) .lt. eMin .or. 
     2        dmgFiberT(k) .eq. one ) then
          statusMp(k) = zero
          endif
c
      nDmg = nDmg + lDmg
c
      endif
c     print*, "rft = ", rft 
cc      print*, "rfc = ", rfc
c      print*, "DmgFiberT = ", DmgFiberT(k)
c      print*, "DmgFiberC = ", DmgFiberC(k)
c      print*, "enomMin = ", enomMin
c      print*, "enomMax = ", enomMax
c      print*, "statusMp(k) = ", statusMp(k)
c      print*,"diseqft = ", diseqft
c       print*,"streqft = ", streqft
c       print*,"stressinatft = ", stressinatft
c       print*,"disinatft = ", disinatft
c       print*,"disfailft = ", disfailft       
c       print*,"diseqfc = ", diseqfc
c       print*,"streqfc = ", streqfc
c       print*,"stressinatfc = ", stressinatfc
c       print*,"disinatfc = ", disinatfc
c       print*,"disfailfc = ", disfailfc
c  
c      print*, "rmt = ", rmt
c      print*, "rmc = ", rmc
c      print*, "dmgMatrixT = ", dmgMatrixT(k)
c      print*, "dmgMatrixC = ", dmgMatrixC(k)
      
c      print*,"diseqmt = ", diseqmt
c       print*,"streqmt = ", streqmt
c       print*,"stressinatmt = ", stressinatmt
c       print*,"disinatmt = ", disinatmt
c       print*,"disfailmt = ", disfailmt
c      print*,"CharLength = ", CharLength(k)
c       print*,"diseqmc = ", diseqmc
c       print*,"streqmc = ", streqmc
c       print*,"stressinatmc = ", stressinatmc
c       print*,"disinatmc = ", disinatmc
c       print*,"disfailmc = ", disfailmc
c
c       print*,"s11 = ", stress(k,i_s33_Xx)
c       print*,"s22 = ", stress(k,i_s33_Yy)
c       print*,"s33 = ", stress(k,i_s33_Zz)
c       print*,"s12 = ", stress(k,i_s33_Xy)
c       print*,"s23 = ", stress(k,i_s33_Yz)
c       print*,"s13 = ", stress(k,i_s33_Zx)
c       print*,"STRAIN11 = ",  STRAIN(K,I_S33_XX)
c       print*,"STRAIN22 = ",  STRAIN(K,I_S33_YY)
c       print*,"STRAIN33 = ",  STRAIN(K,I_S33_ZZ)
c       print*,"STRAIN12 = ",  STRAIN(K,I_S33_XY)
c       print*,"STRAIN23 = ",  STRAIN(K,I_S33_YZ),k, statusMp(k)  
c	 
c       print*,"**********************************"
c
      enddo
c
      return
      end
c
************************************************************ 
*   betaDamping: Add beta damping                          * 
************************************************************ 
      subroutine betaDamping3d ( nblock, 
     *     beta, dtArray, strainInc, sigOld, sigNew, 
     *     statusMp, sigDampOld, sigDampNew ) 
* 
      include 'vaba_param.inc' 
* 
      parameter( 
     *     i_s33_Xx = 1, 
     *     i_s33_Yy = 2, 
     *     i_s33_Zz = 3, 
     *     i_s33_Xy = 4, 
     *     i_s33_Yz = 5, 
     *     i_s33_Zx = 6, 
     *     n_s33_Car = 6 ) 
* 
      dimension  sigOld(nblock,n_s33_Car), 
     *     sigNew(nblock,n_s33_Car), 
     *     strainInc(nblock,n_s33_Car), 
     *     statusMp(nblock), 
     *     sigDampOld(nblock,n_s33_Car), 
     *     sigDampNew(nblock,n_s33_Car)       
* 
      parameter ( zero = 0.d0, one = 1.d0, two=2.0d0, 
     *     half = 0.5d0, third = 1.d0/3.d0 ) 
      parameter ( asmall = 1.d-16 ) 
*     
      betaddt =  beta / dtArray 
* 
      do k =1 , nblock 
         sigDampNew(k,i_s33_Xx) = betaddt * statusMp(k) * 
     *        ( sigNew(k,i_s33_Xx) 
     *        - ( sigOld(k,i_s33_Xx) - sigDampOld(k,i_s33_Xx) ) ) 
         sigDampNew(k,i_s33_Yy) = betaddt * statusMp(k) * 
     *        ( sigNew(k,i_s33_Yy) 
     *        - ( sigOld(k,i_s33_Yy) - sigDampOld(k,i_s33_Yy) ) ) 
         sigDampNew(k,i_s33_Zz) = betaddt * statusMp(k) * 
     *        ( sigNew(k,i_s33_Zz) 
     *        - ( sigOld(k,i_s33_Zz) - sigDampOld(k,i_s33_Zz) ) ) 
         sigDampNew(k,i_s33_Xy) = betaddt * statusMp(k) * 
     *        ( sigNew(k,i_s33_Xy) 
     *        - ( sigOld(k,i_s33_Xy) - sigDampOld(k,i_s33_Xy) ) ) 
         sigDampNew(k,i_s33_Yz) = betaddt * statusMp(k) * 
     *        ( sigNew(k,i_s33_Yz) 
     *        - ( sigOld(k,i_s33_Yz) - sigDampOld(k,i_s33_Yz) ) ) 
         sigDampNew(k,i_s33_Zx) = betaddt * statusMp(k) * 
     *        ( sigNew(k,i_s33_Zx) 
     *        - ( sigOld(k,i_s33_Zx) - sigDampOld(k,i_s33_Zx) ) ) 
* 
         sigNew(k,i_s33_Xx) = sigNew(k,i_s33_Xx)+sigDampNew(k,i_s33_Xx) 
         sigNew(k,i_s33_Yy) = sigNew(k,i_s33_Yy)+sigDampNew(k,i_s33_Yy) 
         sigNew(k,i_s33_Zz) = sigNew(k,i_s33_Zz)+sigDampNew(k,i_s33_Zz) 
         sigNew(k,i_s33_Xy) = sigNew(k,i_s33_Xy)+sigDampNew(k,i_s33_Xy) 
         sigNew(k,i_s33_Yz) = sigNew(k,i_s33_Yz)+sigDampNew(k,i_s33_Yz) 
         sigNew(k,i_s33_Zx) = sigNew(k,i_s33_Zx)+sigDampNew(k,i_s33_Zx) 
* 
      end do 
*     
      return 
      end 
************************************************************ 
*   EnergyInternal3d: Compute internal energy for 3d case  * 
************************************************************ 
      subroutine EnergyInternal3d(nblock, sigOld, sigNew , 
     *   strainInc, curDensity, enerInternOld, enerInternNew) 
* 
      include 'vaba_param.inc' 
* 
      parameter( 
     *     i_s33_Xx = 1, 
     *     i_s33_Yy = 2, 
     *     i_s33_Zz = 3, 
     *     i_s33_Xy = 4, 
     *     i_s33_Yz = 5, 
     *     i_s33_Zx = 6, 
     *     n_s33_Car = 6 ) 
* 
      parameter( two = 2.d0, half = .5d0 ) 
* 
      dimension sigOld (nblock,n_s33_Car), sigNew (nblock,n_s33_Car), 
     *     strainInc (nblock,n_s33_Car), curDensity(nblock), 
     *     enerInternOld(nblock), enerInternNew(nblock) 
* 
      do k = 1, nblock 
         stressPower  = half * ( 
     *        ( sigOld(k,i_s33_Xx) + sigNew(k,i_s33_Xx) ) 
     *        * ( strainInc(k,i_s33_Xx) ) 
     *        +       ( sigOld(k,i_s33_Yy) + sigNew(k,i_s33_Yy) ) 
     *        * ( strainInc(k,i_s33_Yy) ) 
     *        +       ( sigOld(k,i_s33_Zz) + sigNew(k,i_s33_Zz) ) 
     *        * ( strainInc(k,i_s33_Zz) ) 
     *        + two * ( sigOld(k,i_s33_Xy) + sigNew(k,i_s33_Xy) ) 
     *        * strainInc(k,i_s33_Xy) 
     *        + two * ( sigOld(k,i_s33_Yz) + sigNew(k,i_s33_Yz) ) 
     *        * strainInc(k,i_s33_Yz) 
     *        + two * ( sigOld(k,i_s33_Zx) + sigNew(k,i_s33_Zx) ) 
     *        * strainInc(k,i_s33_Zx) ) 
*     
         enerInternNew(k) = enerInternOld(k) + stressPower/curDensity(k)
      end do 
*     
      return   
      end   
************************************************************ 
*   CopyR: Copy from one array to another                  * 
************************************************************ 
      subroutine CopyR(nCopy, from, to ) 
* 
      include 'vaba_param.inc' 
* 
      dimension from(nCopy), to(nCopy) 
* 
      do k = 1, nCopy 
         to(k) = from(k) 
      end do 
* 
      return 
      end 
*********************************************************************  
* eig33Anal: Compute eigen values of a 3x3 symmetric matrix analytically * 
*********************************************************************  
      subroutine eig33Anal( nblock, sMat, eigVal ) 
* 
      include 'vaba_param.inc' 
* 
      parameter(i_s33_Xx=1,i_s33_Yy=2,i_s33_Zz=3 ) 
      parameter(i_s33_Xy=4,i_s33_Yz=5,i_s33_Zx=6 ) 
      parameter(i_s33_Yx=i_s33_Xy ) 
      parameter(i_s33_Zy=i_s33_Yz ) 
      parameter(i_s33_Xz=i_s33_Zx,n_s33_Car=6 ) 
* 
      parameter(i_v3d_X=1,i_v3d_Y=2,i_v3d_Z=3 ) 
      parameter(n_v3d_Car=3 ) 
* 
      parameter ( zero = 0.d0, one = 1.d0, two = 2.d0, 
     *     three = 3.d0, half = 0.5d0, third = one / three, 
     *     pi23 = 2.094395102393195d0, 
     *     fuzz = 1.d-8, 
     *     preciz = fuzz * 1.d4 ) 
* 
      dimension eigVal(nblock,n_v3d_Car), sMat(nblock,n_s33_Car) 
* 
      do k = 1, nblock 
        sh  = third*(sMat(k,i_s33_Xx)+sMat(k,i_s33_Yy)+sMat(k,i_s33_Zz))
        s11 = sMat(k,i_s33_Xx) - sh 
        s22 = sMat(k,i_s33_Yy) - sh 
        s33 = sMat(k,i_s33_Zz) - sh 
        s12 = sMat(k,i_s33_Xy) 
        s13 = sMat(k,i_s33_Xz) 
        s23 = sMat(k,i_s33_Yz) 
* 
        fac  = max(abs(s11), abs(s22), abs(s33)) 
        facs = max(abs(s12), abs(s13), abs(s23)) 
        if( facs .lt. (preciz*fac) ) then 
          eigVal(k,i_v3d_X) = sMat(k,i_s33_Xx) 
          eigVal(k,i_v3d_Y) = sMat(k,i_s33_Yy) 
          eigVal(k,i_v3d_Z) = sMat(k,i_s33_Zz) 
        else 
          q = third*((s12**2+s13**2+s23**2)+half*(s11**2+s22**2+s33**2))
          fac = two * sqrt(q) 
          if( fac .gt. fuzz ) then 
            ofac = two/fac 
          else 
            ofac = zero 
          end if 
          s11 = ofac*s11 
          s22 = ofac*s22 
          s33 = ofac*s33 
          s12 = ofac*s12 
          s13 = ofac*s13 
          s23 = ofac*s23 
          r = s12*s13*s23 
     *         + half*(s11*s22*s33-s11*s23**2-s22*s13**2-s33*s12**2) 
          if( r .ge. one-fuzz ) then 
            cos1 = -half 
            cos2 = -half 
            cos3 = one 
          else if( r .le. fuzz-one ) then 
            cos1 = -one 
            cos2 = half 
            cos3 = half 
          else 
            ang = third * acos(r) 
            cos1 = cos(ang) 
            cos2 = cos(ang+pi23) 
            cos3 =-cos1-cos2 
          end if 
          eigVal(k,i_v3d_X) = sh + fac*cos1 
          eigVal(k,i_v3d_Y) = sh + fac*cos2 
          eigVal(k,i_v3d_Z) = sh + fac*cos3 
        end if 
      end do 
* 
      return 
      end
**************************************************************************************
******
**************************************************************************************
******
************************************Operateur ZZMACCAULEY**************************************
**************************************************************************************
******
      Real*4 FUNCTION ZZMACCAULEY(X)
      Real*4  X
c      print*,"before zz =",X
c
      if (X.lt. 0.D0 ) then
      ZZMACCAULEY = 0.D0
      Else
      ZZMACCAULEY = X
      end if
c      print*,"after zz= ",ZZMACCAULEY
      return
      end
C 
