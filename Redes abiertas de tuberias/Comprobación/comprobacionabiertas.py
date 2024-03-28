import math

###############################################################################
################ COMPROBACIÓN DE DISEÑO PARA REDES ABIERTAS ###################
###############################################################################

##########################################
####### Datos globales del sistema #######
##########################################

rho = 999.1 # Densidad del fluido en kg/m^3
mu = 0.00719  #Viscosidad dinámica 
#nu = mu/rho # Viscosidad cinemática en m^2/s
nu = 1.14e-06 # Viscosidad cinemática en m^2/s


#Tanques

RE1 = [43,0] #Principal que suministra, cota de elevación
RE2 = [34,0]
RE3 = [30,0]
RE4 = [12,0]
RE5 = [0,0]

#Nodos
#Nn = [Z (altura) del nodo, Q (Caudal) demandado en el nodo]

N1 = [35, 0.12]
N2 = [20, 0.1]

#Tuberías

#tub_n=[Start (empieza), End (termina),  L, D, epsilon, K, E (1) o S (-1), No. del nodo con el ]

#Ejemplo 1

tub1=[RE1, N1, 800, 0.4572, 3e-05, 0, 1]
tub2=[N1, RE2, 600, 0.1524, 3e-05, 0, -1]
tub3=[N1, RE3, 700, 0.1524, 3e-05, 0, -1]
tub4=[N1, N2, 550, 0.3048, 3e-05, 0, -1] #Defino los signos para el primer nodo
tub5=[N2, RE4, 300, 0.1524, 3e-05, 0, -1]
tub6=[N2, RE5, 200, 0.1524, 3e-05, 0, -1]


sist1={'tub1':tub1,'tub2':tub2,'tub3':tub3,'tub4':tub4}
sist2={'tub4':tub4, 'tub5':tub5, 'tub6':tub6}

#gravedad m/s^2
g=9.81

# Asumir un valor inicial de f
f = 0.01

# Calcular f utilizando la iteración
error = 1e-6  # Tolerancia de error

iteraQ1=[] #Almacena los caudales que interactúan con el nodo1
iteraZ1=[] #Almacena la diferencia de cotas de las tuberías conectadas al nodo1

iteraQ2=[] #Almacena los caudales que interactúan con el nodo1
iteraZ2=[] #Almacena la diferencia de cotas de las tuberías conectadas al nodo1

i=1
while True:
    tub4[6]=-1
    for tub in sist1.values():

        Ein= tub[0][0]
        Eout= tub[1][0]

        V=math.sqrt((2*g*(Ein-Eout))*(((f*(tub[2]/tub[3]))+tub[5])**(-1.)))
        A = ((math.pi)/4)*(tub[3]**2)

        Re = (V * tub[3]) / nu     
        
        while True:

            t1=tub[4]/(3.7*tub[3])
            t2=2.51/(Re*math.sqrt(f))

            #f_new = (-2 * math.log10((epsilon / (3.7 * D)) + (2.52 / (Re * math.sqrt(f))))) ** (-2)   

            f_new=(-2.*math.log10(t1+t2))**(-2.)

            V_new=math.sqrt((2*g*(Ein-Eout))*(((f_new*(tub[2]/tub[3]))+tub[5])**(-1.)))

            Re_new = (V_new * tub[3]) / nu


            # Verificar la convergencia
            if abs(f - f_new) <= error:
                break 

            f = f_new
            V = V_new
            Re = Re_new

        Q=V*A

        Q=Q*tub[6]

        print(f"Itera {i:.6f}. El caudal (Q) para la tub es: {Q:.6f}")

        iteraQ1.append(Q)

        Z1=Ein-Eout

        print(f"Itera {i:.6f}. La dif de alturas (Z2) para la tub es: {Z1:.6f}")

        iteraZ1.append(Z1)

    #denoC1=sum(iteraQ1/iteraZ1)
    vabs_iteraQ1=[abs(x) for x in iteraQ1]
    a = [x / y for x, y in zip(vabs_iteraQ1, iteraZ1)]
    denoC1=sum(a)
    print(sum(iteraQ1))

    if abs(sum(iteraQ1)-N1[1])<=error:
        break

    else:
        C1=(2*(sum(iteraQ1)-N1[1]))/(denoC1)
        N1[0]=N1[0]+C1

    tub4[6]=1

    for tub in sist2.values():

        Ein= tub[0][0]
        Eout= tub[1][0]

        V=math.sqrt((2*g*(Ein-Eout))*(((f*(tub[2]/tub[3]))+tub[5])**(-1.)))
        A = ((math.pi)/4)*(tub[3]**2)

        Re = (V * tub[3]) / nu     
        
        while True:

            t1=tub[4]/(3.7*tub[3])
            t2=2.51/(Re*math.sqrt(f))

            #f_new = (-2 * math.log10((epsilon / (3.7 * D)) + (2.52 / (Re * math.sqrt(f))))) ** (-2)   

            f_new=(-2.*math.log10(t1+t2))**(-2.)

            V_new=math.sqrt((2*g*(Ein-Eout))*(((f_new*(tub[2]/tub[3]))+tub[5])**(-1.)))

            Re_new = (V_new * tub[3]) / nu


            # Verificar la convergencia
            if abs(f - f_new) <= error:
                break 

            f = f_new
            V = V_new
            Re = Re_new

        Q=V*A

        Q=Q*tub[6]

        print(f"Itera {i:.6f}. El caudal (Q) para la tub es: {Q:.6f}")

        iteraQ2.append(Q)

        Z2=Ein-Eout

        print(f"Itera {i:.6f}. La dif de alturas (Z2) para la tub es: {Z2:.6f}")

        iteraZ2.append(Z2)

    #denoC1=sum(iteraQ1/iteraZ1)
    vabs_iteraQ2=[abs(x) for x in iteraQ2]
    a = [x / y for x, y in zip(vabs_iteraQ2, iteraZ2)]
    denoC2=sum(a)
    print(sum(iteraQ2))

    if abs(sum(iteraQ2)-N2[1])<=error:
        break

    else:
        C2=(2*(sum(iteraQ2)-N2[1]))/(denoC2)
        N2[0]=N2[0]+C2

    iteraQ1=[]
    iteraQ2=[]
    iteraZ1=[]
    iteraZ2=[]
    i=i+1


        


















