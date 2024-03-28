import math

##################################################################################
################ CALCULO DE POTENCIA PARA TUBERÍAS EN PARALELO ###################
##################################################################################

##########################################
####### Datos globales del sistema #######
##########################################


rho = 860.0 # Densidad del fluido en kg/m^3
mu = 0.00719  #Viscosidad dinámica 
nu = mu/rho # Viscosidad cinemática en m^2/s
#nu = 1.17e-06 # Viscosidad cinemática en m^2/s
QT = 0.46 #Caudal total del sistema

E1 = 103.7148 #Energia en el nodo de entrada (kPa)
E2 = 0 #En el nodo se salida


#######################################################
###### Datos locales de cada una de las tuberias ######
#######################################################

#tub_n=[hb/ht, efi, L, D, epsilon, K]

#Ejemplo 1

tub1=[0, 0, 278, 0.4572, 4.6e-05, 7.7]
tub2=[0, 0, 312, 0.3048, 4.6e-05, 9.4]
tub3=[0, 0, 312, 0.3048, 4.6e-05, 9.4]


#Defino la cantidad de tuberías totales en el sistema
sist={'tub1':tub1,'tub2':tub2, 'tub3':tub3}

#gravedad m/s^2
g=9.81

# Asumir un valor inicial de f
f = 0.01

# Calcular f utilizando la iteración
error = 1e-6  # Tolerancia de error


#### Energia
#E1 = (tub1[3]*tub1[4])/(rho*g*Q) #Energia en 1 (Verificar si hay una bomba) o la cabeza 
#E2 = 0 #Energia en 2

iteraQ=[]

j=0
p=0
while True:
    #Q= -1*(tub1[0]+tub2[0]+tub3[0])#Caudal total
    i=0
    for tub in sist.values():

        if j !=0 and i==0:
            Q1=Q1*(QT/QT_new)
            iteraQ=[]

            V=Q1/(((math.pi)/4)*(tub[3]**2)) #Calculo de la velocidad

            Re = (V * tub[3]) / nu

            while True:
                    
                t1=tub[4]/(3.7*tub[3])
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

            hf= f*(tub[2]/tub[3])*(V**2)/(2*g)
            he= tub[5]*((V**2)/(2*g))

            HT=hf+he

            E2=E1-HT

        
        if tub == sist['tub1'] and p==0:
            denoQ1=0
            a=0
            for tuba in sist.values(): #tuba solo itera para obtener la sumatoria de cada una de las longitudes y diametros de las tuberias
                a=(tuba[3]**(5/2))*((math.sqrt(tuba[2]))**(-1.))
                denoQ1+=a #sumatoria

            Q1=QT*(((tub[3]**(5/2))*((math.sqrt(tub[2]))**(-1.)))*((denoQ1)**(-1.)))

            iteraQ.append(Q1)

        #Calculo de la cabeza 

        
            V=Q1/(((math.pi)/4)*(tub[3]**2)) #Calculo de la velocidad

            Re = (V * tub[3]) / nu

            while True:
                    
                t1=tub[4]/(3.7*tub[3])
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

            hf= f*(tub[2]/tub[3])*(V**2)/(2*g)
            he= tub[5]*((V**2)/(2*g))

            HT=hf+he

            E2=E1-HT

        else:
            
            #E2=E1-HT

            V=math.sqrt((2*g*(E1-E2+tub[0]+tub[0]))*(((f*(tub[2]/tub[3]))+tub[5])**(-1.)))
            #V=((-2*math.sqrt(2*g*tub[3]*hf))/(math.sqrt(tub[2])))*(math.log10(((tub[4])/(3.7*tub[3]))+((2.52*nu*math.sqrt(tub[2]))/(tub[3]*math.sqrt(2*g*tub[3]*hf)))))

            Re = (V * tub[3]) / nu

            f = 0.01

            while True:

                t1=tub[4]/(3.7*tub[3])
                t2=2.51/(Re*math.sqrt(f))

                #f_new = (-2 * math.log10((epsilon / (3.7 * D)) + (2.52 / (Re * math.sqrt(f))))) ** (-2)   

                f_new=(-2.*math.log10(t1+t2))**(-2.)

                V_new=math.sqrt((2*g*(E1-E2+tub[0]+tub[0]))*(((f_new*(tub[2]/tub[3]))+tub[5])**(-1.)))

                Re_new = (V_new * tub[3]) / nu

                # Verificar la convergencia
                if abs(f - f_new) <= error:
                    break 

                f = f_new
                V = V_new
                Re = Re_new

            Q=V*(((math.pi)/4)*(tub[3]**2.))

            print(Q)
    
            iteraQ.append(Q)

        i=i+1
        p=p+1
    
    QT_new=sum(iteraQ)

    j=j+1 #contador de iteraciones globales


    if abs(QT-QT_new)<=error: #Convergencia GLOBAL de las iteraciones
        break

P2=rho*g*E2

print(f"La cabeza de presión en el nodo 2 es: {E2:.6f} (m)")

print(f"La presión en el nodo 2 es: {P2:.6f} (Pa)")


    


    




