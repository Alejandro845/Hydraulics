import math

################################################################################################
################ CALCULO DE COMPROBACIÓN DE DISEÑO PARA TUBERÍAS EN PARALELO ###################
################################################################################################

##########################################
####### Datos globales del sistema #######
##########################################


rho = 999.1 # Densidad del fluido en kg/m^3
mu = 0.00114  #Viscosidad dinámica 
nu = mu/rho # Viscosidad cinemática en m^2/s
#nu = 1.14e-06 # Viscosidad cinemática en m^2/s

E1 = 54.0265 #Energia en el sistema en presión normalmente
E2 = 13.2518 #(kPa)

#######################################################
###### Datos locales de cada una de las tuberias ######
#######################################################

#tub_n=[Qo (Caudal en la salida o lateral), Qi (Caudal de entrada de la tuberia), hb o ht (cabeza de bomba o turb), efi (eficiencia de la bomba o turbina), D (diametro), L (long), epsilon (rugosidad), K (sumatoria de perdidas menores)]
#tub_n=[hb/ht, efi, L, D, epsilon, K]

#Ejemplo 1

tub1=[0, 0, 278, 0.1524, 4.6e-05, 7.4]
tub2=[0, 0, 230, 0.2031, 4.6e-05, 6.1]
tub3=[0, 0, 278, 0.1524, 4.6e-05, 7.4]

sist={'tub1':tub1,'tub2':tub2,'tub3':tub3}

#gravedad m/s^2
g=9.81

# Asumir un valor inicial de f
f = 0.01

# Calcular f utilizando la iteración
error = 1e-6  # Tolerancia de error

Q_itera=[]
j=0
delta=0
i=0

for tub in sist.values():

    V=math.sqrt((2*g*(E1-E2+tub[1]+tub[1]))*((f*(tub[2]/tub[3])+tub[5])**(-1.)))

    Re = (V * tub[3]) / nu
        
    while True:
                
        t1=tub[4]/(3.7*tub[3])
        t2=2.51/(Re*math.sqrt(f))

        f_new = (-2.*math.log10(t1+t2))**(-2.)

        V_new = math.sqrt((2*g*(E1-E2+tub[1]+tub[1]))*((f_new*(tub[2]/tub[3])+tub[5])**(-1.)))

        Re_new = (V_new * tub[3]) / nu


        # Verificar la convergencia - condición 1
        if abs(f - f_new) <= error:
            break 

        f = f_new
        V = V_new
        Re = Re_new
 
    print(f"el factor de fric (f) en la itera es: {f_new:.6f}")
    print(f"La velocidad (Vc) en la itera es: {V:.6f}")

    #Calcula el caudal nuevo con la ultima velocidad de la iteracion anterior
    Q=V*((math.pi)/4)*(tub[3]**2.)

    print(f"El caudal para la tubería es: {Q:.6f}")

    Q_itera.append(Q)

    i=i+1

QT=sum(Q_itera)

print(f"El caudal total es: {QT:.6f}")

        