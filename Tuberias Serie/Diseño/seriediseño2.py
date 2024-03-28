import math

#############################################################################
################ CALCULO DE DISEÑO PARA TUBERÍAS EN SERIE ###################
#############################################################################

##########################################
####### Datos globales del sistema #######
##########################################


rho = 999.1 # Densidad del fluido en kg/m^3
#mu = 0.001005  #Viscosidad dinámica 
#nu = mu/rho # Viscosidad cinemática en m^2/s
nu = 1.14e-06 # Viscosidad cinemática en m^2/s
ki = 0 #Coeficiente de perdidas menores en la entrada
ko = 0 #Coeficiente de perdidas menores en la salida


#######################################################
###### Datos locales de cada una de las tuberias ######
#######################################################

#tub_n=[Qo (Caudal en la salida o lateral), Qi (Caudal de entrada de la tuberia), hb o ht (cabeza de bomba o turb), efi (eficiencia de la bomba o turbina), D (diametro), L (long), epsilon (rugosidad), K (sumatoria de perdidas menores)]
#tub_n=[Qo, Qi, hb/ht, efi, P, L, epsilon, K]

#Ejemplo 2

tub1=[-0.16, 0, 0, 0, 0, 72, 0.00015, 6.4]
tub2=[-0.16, 0, 0, 0, 0, 45, 0.00015, 3.6]
tub3=[-0.16, 0, 0, 0, 0, 45, 0.00015, 3.6]



#Diametros comerciales en pulgadas/mm

#CD={'3/4':23.63,'1':30.20,'1 1/4':38.14,'1 1/2':43.68,'2':54.58,'2 1/12':66.07,'3':80.42,'4':103.42,'6':152.22,'8':198.21,'10':247.09,'12':293.07,'14':321.76,'16':367.70,'18':413.66,'20':459.64,'24':551.54}
CD={'2':50.8,'2 1/12':63.5,'3':76.2,'4':101.6,'6':152.4,'8':203.2,'10':254,'12':304.8,'14':355.6,'16':406.4} 
#Creo una lista para almacenar solo los valores de los diametros
diametros=[]

#Ciclo para extraer y almcanear los valores en la lista de diametros
for diam in CD:
    diametros.append(float(CD[diam]))

#Defino la cantidad de tuberías totales en el sistema
sist={'tub1':tub1,'tub2':tub2, 'tub3':tub3}

#gravedad m/s^2
g=9.81

# Asumir un valor inicial de f
f = 0.01

# Calcular f utilizando la iteración
error = 1e-6  # Tolerancia de error

Q= -1*(tub1[0]+tub2[0]+tub3[0])#Caudal total

#### Energia
E1 = 19.3 #Energia en 1 (Verificar si hay una bomba) o la cabeza 
E2 = 0 #Energia en 2

Z1=19.3
Z2=0

iteraL=[] #almacena los valores de longitud de cada tubería

for tub in sist.values(): #pequeña iteración para almacenar las longitudes
    iteraL.append(tub[5])

cosT=math.cos(math.asin((abs(Z1-Z2))/(sum(iteraL)))) #Halla el valor del angulo para con la dif de cabezas y la longitud total

HT=E1-E2

iterahf=[] #Inicializo las listas donde voy a almacenar los datos de las iteraciones
iterahe=[]

hmvi=0 #variable para controlar la convergencia de la energía

j=0
while True:
    Q= -1*(tub1[0]+tub2[0]+tub3[0])#Caudal total
    i=0
    p=0
    for tub in sist.values():
        
        if j == 0: #Calculo del hfi solo pra la primera iteración GLOBAL
            denohf=0
            a=0
            for tuba in sist.values(): #tuba solo itera para obtener la sumatoria de cada una de las longitudes y diametros de las tuberias
                a=(tuba[5])*(cosT)
                denohf+=a #sumatoria

            hfi=HT*((tub[5]*cosT)/(denohf)) #Calculo el hfi para cada tub
        else:
            hfi=u[p]+hmvf

        if tub == sist['tub1']:
            Ein=E1  #La energía de entrada SOLO de la primera tubería será E1
        else: #Sino, será la salida de la tubería anterior
            Ein=Eout 

        if tub == sist['tub1']: #Condicion para la resta de los caudales laterales
            Q=Q
        else:
            l=str(i)
            x='tub'+l
            Q=Q+sist[x][0]
        
        Eout=Ein-hfi

        #A = ((math.pi)/4)*(tub[3]**2)
        for c in range(len(diametros)):
            D=diametros[c]

            D=D/1000

            V=math.sqrt((2*g*(Ein-Eout+tub[2]+tub[2]))*((f*(tub[5]/D)+tub[7])**(-1.)))

            Re = (V * D) / nu

            #While para calcular el valor de punto fijo con punto fijo
            while True:
                
                t1=tub[6]/(3.7*D)
                t2=2.51/(Re*math.sqrt(f))

                f_new=(-2.*math.log10(t1+t2))**(-2.)

                #f_new = (-2 * math.log10((tub[6]/ (3.7 * D)) + (2.52 / (Re * math.sqrt(f))))) ** (-2)   

                V_new = math.sqrt((2*g*(Ein-Eout+tub[2]+tub[2]))*((f_new*(tub[5]/D)+tub[7])**(-1.)))

                Re_new = (V_new * D) / nu
        
                # Verificar la convergencia - condición 1
                if abs(f - f_new) <= error:
                    break 

                
                f = f_new
                V = V_new
                Re = Re_new

            #Calcula el caudal nuevo con la ultima velocidad de la iteracion anterior
            Qnew=V*((math.pi)/4)*(D**2)

            # Verifica la convergencia - condición 2
            if (Qnew>=Q):
                break

        VR=Q/(((math.pi)/4)*(D**2))

        Re = (VR * D) / nu

        while True:
                
            t1=tub[6]/(3.7*D)
            t2=2.51/(Re*math.sqrt(f))

            #f_new = (-2 * math.log10((epsilon / (3.7 * D)) + (2.52 / (Re * math.sqrt(f))))) ** (-2)   

            f_new=(-2.*math.log10(t1+t2))**(-2.)

            #V_new = math.sqrt((2*g*(Ein-Eout+tub[2]+tub[2]))/((f_new*(tub[5]/D))+tub[7]))

            #Re_new = (V_new * D) / nu

                # Verificar la convergencia
            if abs(f - f_new) <= error:
                break 

            f = f_new
            #V = V_new
            #Re = Re_new

        hf= f*(tub[5]/D)*((VR**2)/(2*g))
        he= tub[7]*((VR**2)/(2*g))

        iterahf.append(hf) #Alamacena los datos de las perdidas por fricción de cada iteración
        iterahe.append(he) #Alamacena los datos de las perdidas menores de cada iteración

        

        i=i+1 #Variable contadora para saber la iteración de la tubería y quitar el caudal lateral
        p=p+1

    hmvf=HT-sum(iterahf)-sum(iterahe)
    
    print(f"El balance de energía al final de la iteración fue: {hmvf:.6f}")

    if abs(hmvi-hmvf)<=error or hmvf<0 or j>10:
        break
    
    hmvi = hmvf

    u = iterahf
    v = iterahe

    iterahf=[]
    iterahe=[]

    j=j+1




