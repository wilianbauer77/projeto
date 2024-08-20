import numpy
import mpmath
import matplotlib.pyplot as plt
#import pylab

#for python 2.7


#######################

#Número de espécies
Nspecies=2
#kappa parameter
k=numpy.zeros(Nspecies)
alpha=numpy.zeros(Nspecies)
#######################


npara=numpy.zeros(Nspecies,dtype='i4')
nperp=numpy.zeros(Nspecies,dtype='i4')

vparamin=numpy.zeros(Nspecies)
vparamax=numpy.zeros(Nspecies)

vperpmin=numpy.zeros(Nspecies)
vperpmax=numpy.zeros(Nspecies)

dens=numpy.zeros(Nspecies)
mu=numpy.zeros(Nspecies)

beta_para=numpy.zeros(Nspecies)
beta_perp=numpy.zeros(Nspecies)

theta_para=numpy.zeros(Nspecies)
theta_perp=numpy.zeros(Nspecies)
vdrift=numpy.zeros(Nspecies)

########################

#species 1
k[0]=4
alpha[0]=0.05

npara[0]=255
nperp[0]=64

dens[0]=1.0
mu[0]=1.0

beta_para[0]=4.0
beta_perp[0]=2.0

theta_para[0]=(((k[0]-1.5)/k[0])*beta_para[0]*mu[0])**(0.5)
theta_perp[0]=(((k[0]-1.5)/k[0])*beta_perp[0]*mu[0])**(0.5)

vdrift[0]=-0.0

vparamin[0]=-12.0
vparamax[0]=12.0

vperpmin[0]=0.0
vperpmax[0]=10.0


#species 2
k[1]=50
alpha[0]=0.2

npara[1]=255
nperp[1]=64

dens[1]=1.0
mu[1]=1836.0

beta_para[1]=1.0
beta_perp[1]=1.0

theta_para[1]=(((k[1]-1.5)/k[1])*beta_para[1]*mu[1])**(0.5)
theta_perp[1]=(((k[1]-1.5)/k[1])*beta_perp[1]*mu[1])**(0.5)

vdrift[1]=-0.0

vparamin[1]=-260.0
vparamax[1]=260.0

vperpmin[1]=0.0
vperpmax[1]=260.0

#############################

limit=10.0**(-300)


def dist_bikappa(vpar, vper, n,m, kappa,a,theta_par,theta_per,drift):
	theta_par_per2 = numpy.outer(theta_par, theta_per**2)
	theta3=theta_par*theta_per**2
	N=(1/((numpy.pi*kappa)**(1.5)*theta3))*mpmath.hyp1f1(1.5,1.5-kappa,a*kappa)
	bikappa=(n/(kappa**1.5*numpy.pi**1.5*theta_par_per2))*(mpmath.gamma(kappa+1)/(mpmath.gamma(kappa-0.5)))*(1+vpar**2/(kappa*theta_par**2)+vper**2/(kappa*theta_per**2))**(-kappa-1)
	#return bikappa
	return bikappa.ravel()

for ispecies in range(0,Nspecies):

	file_name='distribution_bikappa'+str(ispecies+1)+'.dat'

	vpara = numpy.linspace(vparamin[ispecies],vparamax[ispecies], npara[ispecies])
	vperp = numpy.linspace(vperpmin[ispecies],vperpmax[ispecies], nperp[ispecies])

	vpara2,vperp2=numpy.meshgrid(vpara,vperp)

	data=dist_bikappa(vpara2, vperp2,dens[ispecies],mu[ispecies],k[ispecies],theta_para[ispecies],theta_perp[ispecies],vdrift[ispecies])
	data=data.reshape(nperp[ispecies],npara[ispecies])

	data_new=numpy.zeros((nperp[ispecies],npara[ispecies]))

	for i in range(0,npara[ispecies]):
		for j in range(0,nperp[ispecies]):

			if(data[j,i]>limit):
				data_new[j,i]=data[j,i]*dens[ispecies]
			else:
				data_new[j,i]=0.0


	dat_fin = open(file_name, 'w')

	for i in range(0,npara[ispecies]):
		for j in range(0,nperp[ispecies]):
			dat_fin.write(str(vpara[i]))
			dat_fin.write(" ")
			dat_fin.write(str(vperp[j]))
			dat_fin.write(" ")
			dat_fin.write(str(data_new[j,i]))
			dat_fin.write("\n")
	
	
	# Plotagem do gráfico
	plt.plot(vpara, data_new[:, :].sum(axis=0))  # Soma sobre a dimensão vperp
	plt.xlabel('vpara')
	plt.ylabel('fk(v)')
	plt.title('Species ' + str(ispecies + 1) + ' Distribution')
	plt.savefig(file_name + '.png')  # Salva o gráfico como imagem

	plt.show()  # Exibe todos os gráficos
