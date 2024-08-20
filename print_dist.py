import numpy
import mpmath
import matplotlib.pyplot as plt
import array as arr
#import pylab

#for python 2.7

#######################

#Número de espécies
Nspecies=2
#kappa parameter
nkappa=5
k=numpy.zeros(Nspecies)
#K=[Nspecies][nkappa]
K=[[0]*nkappa for i in range(Nspecies)]

K = [[2,4,10,100,1000],[2,4,10,100,1000]]




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
k=K[0][:]

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

k=K[1][:]

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

def dist_bimax(vpar, vper, n,m,beta_par,beta_per,drift):
	bimax=numpy.exp(-n*(vpar-drift)**2/beta_par/m -n*vper**2/beta_per/m)* n**1.5 /(m**1.5 *numpy.pi**1.5 *beta_per*numpy.sqrt(beta_par))
	return bimax.ravel()

def dist_bikappa(vpar, vper, n, kappa,theta_par,theta_per,drift):
	theta_par_per2 = numpy.outer(theta_par, theta_per**2)
	bikappa=(n*(1.5)/(kappa**1.5*numpy.pi**1.5*theta_par_per2))*(mpmath.gamma(kappa+1)/(mpmath.gamma(kappa-0.5)))*(1+(n)*(vpar-drift)**2/(kappa*theta_par**2)+(n)*vper**2/(kappa*theta_per**2))**(-kappa-1)
	#return bikappa
	return bikappa.ravel()

for ispecies in range(0,Nspecies):
	file_name='distribution'+str(ispecies+1)+'.dat'
	#file_name_bimax='distribution_bimax'+str(ispecies+1)+'.dat'
	#file_name_bikappa='distribution_bikappa_'+str(ispecies+1)+'_'+str(K[ispecies][ki])+'.dat'

	vpara = numpy.linspace(vparamin[ispecies],vparamax[ispecies], npara[ispecies])
	vperp = numpy.linspace(vperpmin[ispecies],vperpmax[ispecies], nperp[ispecies])

	vpara2,vperp2=numpy.meshgrid(vpara,vperp)

	databimax=dist_bimax(vpara2, vperp2,dens[ispecies],mu[ispecies],beta_para[ispecies],beta_perp[ispecies],vdrift[ispecies])
	databimax=databimax.reshape(nperp[ispecies],npara[ispecies])
		

	#data_new=numpy.zeros((nperp[ispecies],npara[ispecies]))
	data_new_bimax=numpy.zeros((nperp[ispecies],npara[ispecies]))
	data_new_bikappa=numpy.zeros((nperp[ispecies],npara[ispecies]))
	for i in range(0,npara[ispecies]):
		for j in range(0,nperp[ispecies]):
			if(databimax[j,i] >limit):
				data_new_bimax[j,i]=databimax[j,i]*dens[ispecies]
			else:
				data_new_bimax[j,i]=0.0

	
	for ki in range(0,nkappa):
		databikappa=dist_bikappa(vpara2, vperp2,dens[ispecies],K[ispecies][ki],theta_para[ispecies],theta_perp[ispecies],vdrift[ispecies])
		databikappa=databikappa.reshape(nperp[ispecies],npara[ispecies])
		for i in range(0,npara[ispecies]):
			for j in range(0,nperp[ispecies]):
				if(databikappa[j,i] >limit):
					data_new_bikappa[j,i]=databikappa[j,i]*dens[ispecies]
				else:
					data_new_bikappa[j,i]=0.0


		#dat_fin_bimax = open(file_name_bimax, 'w')
		#dat_fin_bikappa = open(file_name_bikappa, 'w')
	dat_fin = open(file_name, 'w')


	for i in range(0,npara[ispecies]):
		for j in range(0,nperp[ispecies]):
			dat_fin.write(str(vpara[i]))
			dat_fin.write(" ")
			dat_fin.write(str(vperp[j]))
			dat_fin.write(" ")
			dat_fin.write(str(data_new_bimax[j,i]))
			dat_fin.write(" ")
			dat_fin.write(str(data_new_bikappa[j,i]))
			dat_fin.write("\n")
		
		
		# Plotagem do gráfico
		plt.plot(vpara, numpy.log(data_new_bimax[:, :]).sum(axis=0), color='r', label='maxwellian')
		plt.plot(vpara, numpy.log(data_new_bikappa[:, :]).sum(axis=0),color='g', label='kappa='+str(K[ispecies][ki]))
		# Soma sobre a dimensão vperp
		plt.xlabel('vpara')
		plt.ylabel('f(v)')
		plt.title('Species ' + str(ispecies + 1) + ' Distribution')
		plt.legend()
		#plt.savefig(file_name_bimax + '.png')  # Salva o gráfico como imagem
		#plt.savefig(file_name_bikappa + '.png')  # Salva o gráfico como imagem
		plt.savefig(file_name + '.png')
		plt.show()  # Exibe todos os gráficos
