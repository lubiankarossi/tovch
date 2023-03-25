
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!																			   !!!
!!						TOV Charged Equations Solver						   !!!
!!																			   !!!
!!					   Autores:												   !!!
!!					           P.C.R.Rossí									   !!!
!!					           L.F.R.Rossi									   !!!
!!																			   !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!	  																		   !!
!!	Este programa resolve as equações de equilibrio hidrostático TOV para es - !!
!!	trelas carregadas, usando o método de Runge-Kutta de quinta ordem. 		   !!
!!																			   !!
!!	As grandezas fisicas estão em unidades naturais (c=1,G=1).				   !!
!!																			   !!
!!	Indice de variaveis :													   !!
!!																			   !!
!!				t    = raio (r)												   !!
!!				y(1) = pressão (P)											   !!
!!				y(2) = curvatura polar (lambda) 							   !!
!!				y(3) = carga (Q)											   !!
!!				y(4) = massa (M)											   !!
!!				y(5) = curvatura temporal (v)								   !!
!!				y(6) = densidade de energia (e)								   !!
!!																			   !!
!!						 ou													   !!
!!																			   !!
!!				Y = [P,Lambda,Q,M,v,e]										   !!
!!																			   !!
!!	O Programa chama quatro arquivos esternos de entrada("inp.dat","scale.dat",!!
!!  "config.dat","EOS.dat").															   !!		   
!!																			   !!
!! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!											   !!
!! !!						   !!											   !!
!! !!	Descrição de "inp.dat" !!											   !!
!! !!						   !!											   !!
!! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!											   !!
!!																			   !!
!!	Este arquivo é uma listagem de 4 dados em linhas diferentes, os primeiros  !!
!!	dois dados "pin e pf" especificam uma faixa de pressão [pin,pf] para o cal-!!
!!	culo (pin<pf), o terceiro "n" indica o numero de pontos considerado nesta  !!
!!  faixa e igualmente distribuidos, n é um inteiro maior ou iqual a 1, se n=1 !!
!!	então pin é considerado e pf descartado, o quarto "flag" aceita somente 2  !!
!!	números inteiros, ou 0 ou 1, se flag=0 obtem a densidade de energia da fun-!! 
!!	ção de estado arbitraria e=e(t,y) a ser especificada logo abaixo (no começo!! 
!!	do bloco de subrotinas e funções e se flag=1 obtem a densidade de energia  !!
!!	em um arquivo de dados "EOS.DAT" cuja a forma esta descrita abaixo.		   !!
!!																			   !!
!!			 forma do inp.dat :												   !!
!!																			   !!
!!				pin  - pressão central inferior								   !!
!!				pf   - pressão central superior								   !!
!!				n    - numero de pontos considerado em [pin,pf]				   !!
!!				flag - 0 ou 1 (função/dados)								   !!
!!																			   !!
!! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!											   !!
!! !!						   !!											   !!
!! !!	Descrição de "EOS.dat" !!											   !!
!! !!						   !!											   !!
!! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!											   !!
!!																			   !!
!! Este arquivo somente é aplicavel se flag=1 em "inp.dat". Especifica uma lis-!!
!! composta de duas colunas - pressão e densidade de energia, respectivamente -!!
!! em unidades naturais, interpolando ponto a ponto a equação de estado do sis-!!
!! tema. 																	   !!
!!																			   !!
!!					  EOS.DAT {P[fm^-4],erg[fm^-4]}							   !!
!!																			   !!
!! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!											   !!
!! !!						   !!											   !!
!! !! Descrição de "scale.dat" !!											   !!
!! !!						   !!											   !!
!! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!											   !!
!!																			   !!
!! Este arquivo consiste de 7 entrada A,B,....,G que transformam dados de saida!!
!! segundo a lei															   !!
!!																			   !!
!!				t        = 	A * t											   !!
!!				Y(1)	 =	B * Y(1)										   !!
!!				Y(2)	 =	C * Y(2)										   !!
!!				Y(3)	 =	D * Y(3)										   !!
!!				Y(4)	 =	E * Y(4)										   !!
!!				Y(5)	 =	F * Y(5)										   !!
!!				Y(6)	 =	G * Y(6)										   !!
!!																			   !!
!! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!											   !!
!! !!						    !!											   !!
!! !! Descrição de "config.dat" !!											   !!
!! !!						    !!											   !!
!! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!											   !!
!!																			   !!
!! Arquivo de configuração do programa especifica algumas variaveis as quais   !!
!! são : 																	   !!
!!																			   !!
!!				h    - passo inicial do metodo								   !!
!!				dr   - intervalo de amostragem (aplicavel se n=1)			   !!
!!				zero - real que o programa entende por zero					   !!
!!																			   !!
!!																			   !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!																			   !!
!!							SAIDAS DE DADOS									   !!
!!																			   !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!																			   !!
!!	Se n = 1 a saida é o perfil da estrela cuja pressão central é pin, e os ar-!!
!!	quivos de saida são: 			   										   !!
!!																			   !!
!!																			   !!
!!		"es_rpem.out" -> lista -> [r,P,e,M]  								   !!
!!																			   !!
!!			 e																   !!
!!																			   !!
!!		"es_rlvq.out" -> lista -> [r,lambda,v,Q]							   !!
!!																			   !!
!!	 Já se n != 1 então é uma lista de valores extremos para a região [pin,pf] !!
!!	 e cuja saida é															   !!
!!																			   !!
!!		"es_range.out"-> lista -> [r*,pc,M*,Q*]								   !!
!!																			   !!
!!	 onde 																	   !!
!!																			   !!
!!		r* = raio da estrela cuja pressão centel é pc						   !!
!!		pc = pressão central												   !!
!!		M* = massa encerrada no raio r*										   !!
!!		Q* = carga encerrada no raio r*										   !!
!!																			   !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!																			   !!
!!					Definição da densidade de carga "pch"					   !!
!!																			   !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!																			   !!
!!	A densidade de carga deve se definida na função pch(t,y) no inicio do bloco!! 
!!	de subrotinas e funções abaixo.											   !!
!!																			   !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


program tovch_princ
integer, parameter :: dim = 6;
double precision :: y(dim),t;
double precision :: pin,pf,npoints,h,dr,zero;
double precision :: a_,b_,c_,d_,e_,f_,g_;
double precision :: dp,ddr;
integer :: flag,j;

Open(9,File='config.dat',status='unknown') ! abre arquivos de configuração
Open(10,File='inp.dat',status='unknown')   ! abre arquivos de entrada     
Open(11,File='EOS.dat',status='unknown')   ! abre arquivos de dados de estado (EOS)
Open(12,File='scale.dat',status='unknown') ! abre arquivos de escala      

print*,'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
print*,'!!                                                         !!'
print*,'!!                 Runge-Kutta Solver                      !!'
print*,'!!  Solver System of 5 Ordinary Differential Equations     !!'
print*,'!!                  TOV charged Equation                   !!'
print*,'!!                                                         !!'
print*,'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'

print*,'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
print*,'!!                                                         !!'
print*,'!!                      RUNNING !!!                        !!'
print*,'!!                                                         !!'
print*,'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'

Read(9,*)h    ! lê o passo
Read(9,*)dr   ! lê o intervalo de amostragem
Read(9,*)zero ! lê o zero

Read(10,*)pin	   ! lê o pressão inicial
Read(10,*)pf	   ! lê o pressão final
Read(10,*)npoints  ! lê o número de estrelas consideradas
Read(10,*)flag	   ! lê o flag 

Read(12,*)a_  	   ! lê as escalas
Read(12,*)b_  	   !       ||
Read(12,*)c_  	   !       ||
Read(12,*)d_  	   !       ||       
Read(12,*)e_  	   !       ||        
Read(12,*)f_  	   !       ||      
Read(12,*)g_  	   !       ||        

close(9)	 ! fecha config.dat
close(10)	 ! fecha inp.dat
close(12)	 ! fecha scale.dat

dp = (pf-pin)/npoints ! calcula intervalo de progressão em p

if (npoints == 1) then	! se é apenas uma estrela entra
                        !    para calculo de perfil

	Open(2,File='es_rpem.out',status='unknown')	 ! abre\cria arquicos de saidas
	Open(3,File='es_rlvq.out',status='unknown')	 !		   p\ n=1

	t	 = zero			   ! seta valores de partida e tempo
	y(1) = pin			   !	do siatema de equações 
	y(2) = 0.0d0		   !
	y(3) = 0.0d0		   !
	y(4) = 0.0d0		   !
	y(5) = 0.0d0		   !
	y(6) = perg(t,y,flag)  !

	ddr = 0
	print*,'R in [Km] : ' 

	do
		t = t + h
		ddr = ddr + h
		call rk54(y,t,h,flag)

		if (ddr >= dr) then
			write(2,*) Sngl(a_*t),Sngl(b_*y(1)),Sngl(g_*y(6)),Sngl(e_*y(4))
			write(3,*) Sngl(a_*t),Sngl(c_*y(2)),Sngl(f_*y(5)),Sngl(d_*y(3))
			print*,Sngl(a_*t)
			ddr = 0
		endif

		if (y(1) < zero)	exit

	end do
			write(2,*) Sngl(a_*t),Sngl(b_*y(1)),Sngl(g_*y(6)),Sngl(e_*y(4))
			write(3,*) Sngl(a_*t),Sngl(c_*y(2)),Sngl(f_*y(5)),Sngl(d_*y(3))
			print*,Sngl(a_*t) 

	close(2)
	close(3)

else

	Open(4,File='es_range.out',status='unknown')
    
	print*,'	R		P		M'
	do j=1,npoints

		pin = pin + dp

		t	 = zero
		y(1) = pin
		y(2) = 0.0d0
		y(3) = 0.0d0
		y(4) = 0.0d0
		y(5) = 0.0d0
		y(6) = perg(t,y,flag)

		do
			t = t + h
			call rk54(y,t,h,flag)
			if (y(1) < zero) exit
		end do

		write(4,*) Sngl(a_*t),Sngl(b_*pin),Sngl(e_*y(4)),Sngl(d_*y(3))
		print*, Sngl(a_*t),Sngl(b_*pin),Sngl(e_*y(4))

	end do

	close(4)
end if

print*,'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
print*,'!!                                                         !!'
print*,'!!                      FINISHED !!!                       !!'
print*,'!!                                                         !!'
print*,'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'

close(11)



end program tovch_princ

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!																		!!!
!!             Definição das funções desidade de carga e energia		!!!
!!																		!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!																		!!
!!					Função densidade de carga							!!
!!																		!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



function pch(t,y)
integer, parameter :: dim = 6;
double precision :: y(dim),t,alpha;

  alpha   = 0.95;

  pch = alpha*Y(6) 

end function pch


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!																		!!
!!		Função densidade de emergia - aplicavel se flag = 0				!!
!!																		!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

function e(t,y)
integer, parameter :: dim = 6;
double precision :: y(dim),t;

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!															  !!
!					Caso n Relativistico					  !!
!						n=3/2								  !!
!						k=0.05								  !!
!															  !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Double Precision  k_
! Double Precision  gama;
!
!	gama = 5.0d0/3.0;
!	k_   = .05d0;
!
!	e =  abs(y(1)/k_)**(1.0d0/gama)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!															  !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!															  !!
!						a1s0								  !!
!															  !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Double Precision  k_1,k_2,k_3,k_4,k_5,k_
! Double Precision  gama1,gama2,gama3,gama4,gama5,gama;
! Double Precision  p12,p23,p34,p45
! 
! 	k_1 =   exp(-2.9151620d0)        
! 	k_2	=   exp(-2.8987670d0)          
! 	k_3	=   exp(-3.27590d0)        
! 	k_4	=   exp(-1.6733560d0)          
! 	k_5	=   exp(-2.9058440d0)          
!  
! 	gama1 = 3.379436d0        
! 	gama2 = 1.583672d0        
! 	gama3 = 1.931565d0        
! 	gama4 = 1.00d0      
! 	gama5 = 1.624754d0          
! 
! 	p12	=((k_1**gama2)/(k_2**gama1))**(1.0d0/(gama2-gama1))
! 	p23	=((k_2**gama3)/(k_3**gama2))**(1.0d0/(gama3-gama2))
! 	p34	=((k_3**gama4)/(k_4**gama3))**(1.0d0/(gama4-gama3))
! 	p45	=((k_4**gama5)/(k_5**gama4))**(1.0d0/(gama5-gama4))
! 
! 	
! 	if ((y(1) >= 0) .and. (y(1) < p12) ) then 
! 		gama = gama1
! 		k_   = k_1
! 	elseif ((y(1) >= p12) .and. (y(1) < p23)) then
! 		gama = gama2
! 		k_   = k_2
! 	elseif ((y(1) >= p23) .and. (y(1) < p34)) then
! 		gama = gama3
! 		k_   = k_3
! 	elseif ((y(1) >= p34) .and. (y(1) < p45)) then
! 		gama = gama4
! 		k_   = k_4
! 	elseif (y(1) >= p45) then
! 		gama = gama5
! 		k_   = k_5
! 	end if
! 
! 	e =  abs(y(1)/k_)**(1.0d0/gama)
! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!															  !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!															  !!
!						nlwmcfl								  !!
!															  !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!
!Double Precision  k_1,k_2,k_3,k_4,k_
!Double Precision  gama1,gama2,gama3,gama4,gama;
!Double Precision  p12,p23,p34
!
!	k_1 = exp(-2.956244)       
!	k_2	= exp(-2.345218)         
!	k_3 = exp(-5.53365)       
!	k_4	= exp(-2.766383)     
!
!	gama1 = 3.026653      
!	gama2 = 1.024835      
!	gama3 = 2.897337   
!	gama4 =	1.482865
!
!	p12	=(((k_1**gama2)/(k_2**gama1))**(1.0d0/(gama2-gama1)))
!	p23	=(((k_2**gama3)/(k_3**gama2))**(1.0d0/(gama3-gama2)))
!	p34	=(((k_3**gama4)/(k_4**gama3))**(1.0d0/(gama4-gama3)))
!
!	if ((y(1) >= 0) .and. (y(1) < p12) ) then 
!		gama = gama1
!		k_   = k_1
!	elseif ((y(1) >= p12) .and. (y(1) < p23)) then
!		gama = gama2
!		k_   = k_2
!	elseif ((y(1) >= p23) .and. (y(1) < p34)) then
!		gama = gama3
!		k_   = k_3
!	elseif (y(1) >= p34) then
!		gama = gama4
!		k_   = k_4
!	end if
!
!	e = abs(y(1)/k_)**(1/gama)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!															  !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!															  !!
!						pns0								  !!
!															  !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
Double Precision  k_1,k_2,k_3,k_4,k_
Double Precision  gama1,gama2,gama3,gama4,gama;
Double Precision  p12,p23,p34

	k_4 = .3878685028e-1       
	k_3	= .6221868017e-1         
	k_2	= .2122479738       
	k_1	= .7160438241e-1     

	gama4 = 2.3412790      
	gama3 = 1.8547330      
	gama2 = 1.000   
	gama1 =	1.5041720


	p12	=((k_1**gama2)/(k_2**gama1))**(1.0d0/(gama2-gama1))
	p23	=((k_2**gama3)/(k_3**gama2))**(1.0d0/(gama3-gama2))
	p34	=((k_3**gama4)/(k_4**gama3))**(1.0d0/(gama4-gama3))

	
	if ((y(1) >= 0.0) .and. (y(1) < p12) ) then 
 		gama = gama1
		k_   = k_1
	elseif ((y(1) >= p12) .and. (y(1) < p23)) then
		gama = gama2
		k_   = k_2
	elseif ((y(1) >= p23) .and. (y(1) < p34)) then
		gama = gama3
		k_   = k_3
	elseif (y(1) >= p34) then
		gama = gama4
		k_   = k_4
 	end if

	e = abs(y(1)/k_)**(1/gama)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!															  !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
end function e


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!																		!!!
!!	             Rotinas para o metodo e funções   						!!!
!!																		!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!																		!!
!! rk54(y,t,h,flag)  ->  subrotina que efetua a integração pelo método  !!
!!						 de runge-kutta de quinta ordem progredindo em	!!
!!						 t+h os valores de y							!!
!!																		!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine rk54(y,t,h,flag)

integer :: i,flag;
integer, parameter :: dim = 6;
double precision :: y(dim),yt(dim),k(6,dim-1),h,t;
double precision, parameter :: b21 = 0.20d0, &
b31 = 3.0d0/40.0d0,&
b32 = 9.0d0/40.0d0,&
b41 = 0.30d0,&
b42 = -0.90d0,&
b43 = 1.20d0,&
b51 = -11.0d0/54.0d0,&
b52 = 2.50d0,&
b53 = -70.0d0/27.0d0,&
b54 = 35.0d0/27.0d0,&
b61 = 1631.0d0/55296.0d0,&
b62 = 175.0d0/512.0d0,&
b63 = 575.0d0/13824.0d0,&
b64 = 44275.0d0/110592.0d0,&
b65 = 253.0d0/4096.0d0,&
a2  = 0.20d0,&
a3  = 0.30d0,&
a4  = 0.60d0,&
a5  = 1.0d0,&
a6  = 0.875d0,&
c1  = 37.0d0/378.0d0,&
c2  = 0.0d0,&
c3  = 250.0d0/621.0d0,&
c4  = 125.0d0/594.0d0,&
c5  = 0.0d0,&
c6  = 512.0d0/1771.0d0,&
c1_ = c1-2825.0d0/27648.0d0,&
c2_ = 0.0d0,&
c3_ = c3-18575.0d0/48384.0d0,&
c4_ = c4-13525.0d0/55296.0d0,&
c5_ = -277.0d0/14336.0d0,&
c6_ = c6-0.250d0;

yt=y

do i=1,dim-1
	k(1,i) = h * f(t,yt,i)
	yt(i) = y(i) + b21 * k(1,i)
end do
do i=1,dim-1
	k(2,i) = h * f(t+a2*h,yt,i)
	yt(i) = y(i) + b31 * k(1,i) + b32 * k(2,i)
end do
do i=1,dim-1
	k(3,i) = h * f(t+a3*h,yt,i)
	yt(i) = y(i) + b41 * k(1,i) + b42 * k(2,i) + b43 * k(3,i)
end do
do i=1,dim-1
	k(4,i) = h * f(t+a4*h,yt,i)
	yt(i) = y(i) + b51 * k(1,i)+ b52 * k(2,i) + b53 * k(3,i) + b54 * k(4,i)
end do
do i=1,dim-1
	k(5,i) = h * f(t+a5*h,yt,i)
	yt(i) = y(i) + b61 * k(1,i) + b62 * k(2,i) + b63 * k(3,i) + b64 * k(4,i) + b65 * k(5,i)
end do
do i=1,dim-1
	k(6,i) = h * f(t+a6*h,yt,i)
end do
do i=1,dim-1
	y(i) = y(i) + c1*k(1,i)+ c2*k(2,i)+ c3*k(3,i)+ c4*k(4,i)+ c5*k(5,i)+ c6*k(6,i)
end do

y(6) = perg(t+h,y,flag)

end subroutine rk54


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!																		!!
!! f(t,y,index)  -> define as equações diferenciais usadas pelo metodo	!!
!!					onde o indice indexe se refere a cada uma das deri-	!!
!!					vadas												!!
!!																		!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


function f(t,y,index)
integer index;
integer, parameter :: dim = 6;
double precision, parameter :: pi  = 3.1415926535897932
double precision :: y(dim),t,temp,result;

if (index==1) then
	
	result	= -(y(4)+4*pi*(t**3)*(y(1)-(y(3)**2)/(8*pi*(t**4))))*(y(1)+y(6))/ &
	                (t**2-2*y(4)*t) + (y(3)/(t**2))*pch(t,y)*exp(y(2)/2.0)

elseif (index==2) then

	result = 8*pi*t*exp(y(2))*(y(6)+(y(3)**2)/(8*pi*(t**4)))+(1-exp(y(2)))/t

elseif (index==3) then

	result = 4*pi*(t**2)*pch(t,y)*exp(y(2)/2)

elseif (index==4) then

	result	= (4*pi*(t**2)*y(6)+(y(3)**2/(2*t**2)))

elseif (index==5) then

	result	= 8*pi*t*exp(y(2))*(y(1)-(y(3)**2)/(8*pi*(t**4)))+((exp(y(2))-1)/t)

else

	print*, 'Erro em indice de f(t,y,index)'	

endif

  f = result;

end function f


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!																		!!
!!	perg(t,y,flag)  ->  escolhe se a equação de estado é uma função		!!
!!	                    definida em e(t,y) [flag=0]	ou aquivo de		!!
!!						de dados [flag=1]								!!
!!																		!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


function perg(t,y,flag)
integer, parameter :: dim = 6;
double precision :: y(dim),t, result;
integer :: flag;

if (flag == 0) then

	result = e(t,y)	

elseif (flag == 1) then

	call interpol(11,y(1),result)

else

	print*, 'Erro em flag em perg(t,y,flag)'	

endif

perg = result


end function perg


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!																		!!!
!!             Rotinas de Manipulação de Arquivos						!!!
!!																		!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!																		!!
!!	interpol(d,value,result)  ->  interpola valores no arquivo d	    !!
!!																		!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


subroutine interpol(d,X,Y) 
doubleprecision X1,Y1,X2,Y2,a,b,Y,X
integer d,index1,index2


	call findv(d,X,index1,index2)

	if (index1 == index2) then
		
		call get(d,index2,X1,Y1);
		Y = Y1

	else 

		call get(d,index1,X1,Y1)
		call get(d,index2,X2,Y2)
		a =	(Y1-Y2)/(X1-X2);
		b = (Y2*X1-Y1*X2)/(X1-X2);
		Y = a*X+b

	end if 

end subroutine interpol



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!																		!!
!!	findv(x,value,index1,index2)  ->  retorna indices para value no	    !!
!!									  arquivo x na posição index 		!!
!!									  ou faixa (index1,index2)			!!
!!																		!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


subroutine findv(x,value,index1,index2)
doubleprecision count,value
integer x,index1,index2

	REWIND (x)
	count=0;
    DO WHILE (.NOT. EOF(x))
		 read(x,*) value1,value2;
         count = count + 1
		 if (value1 < value) then 
			index1 = count
		 else if (value1.eq.value) then
			index1 = count
			index2 = count
			exit;
		 else 
			index2 = count
			exit;
		 end if

	END DO

end subroutine findv


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!																		!!
!!	get(x,index,value1,value2)	->  retorna valores value1 e 2 do 		!!
!!									  arquivo x na posição index 		!!
!!																		!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


subroutine get(x,index,value1,value2)
doubleprecision count,value1,value2
integer x,index

	REWIND (x)
	count=0;
    DO WHILE (.NOT. EOF(x))
		 read(x,*) value1,value2;
         count = count + 1
    	 if (count.eq.index) exit;
	END DO

end subroutine get


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!																		!!
!!		seek_max(x,count)	->  retorna numero maximo de registos		!!
!!									em count no arquivo x	 			!!
!!																		!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


subroutine seek_max(x,count)
doubleprecision count
integer x

	count=0;
    DO WHILE (.NOT. EOF(x))
		 read(x,*) value1,value2;
         count = count + 1
    
	END DO
end subroutine seek_max

