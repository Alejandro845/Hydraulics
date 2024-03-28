import math

#############################################################################
########## CALCULO DE LA POTENCIA REQUERIDA PARA TUBERÍAS EN SERIE ##########
#############################################################################

##########################################
####### Datos globales del sistema #######
##########################################


rho = 999.3 # Densidad del fluido en kg/m^3
mu = 0.001005  #Viscosidad dinámica 
#nu = mu/rho # Viscosidad cinemática en m^2/s
nu = 1.14e-06 # Viscosidad cinemática en m^2/s
ki = 0 #Coeficiente de perdidas menores en la entrada
ko = 0 #Coeficiente de perdidas menores en la salida
E1 = 0  #Energia en 1
E2 = 31.7 #Energia en 2


#######################################################
###### Datos locales de cada una de las tuberias ######
#######################################################

#tub_n=[Qo (Caudal en la salida o lateral), Qi (Caudal de entrada de la tuberia), hb o ht (cabeza de bomba o turb), efi (eficiencia de la bomba o turbina), D (diametro), L (long), epsilon (rugosidad), K (sumatoria de perdidas menores)]
#tub_n=[Qo, Qi, hb/ht, efi, D, L, epsilon, K]

#Ejemplo 1

tub1=[-0.094, 0, 0, 0.75, 0.2032, 184, 4.6e-05, 7.1]
tub2=[-0.087, 0, 0, 0, 0.1524, 393, 1.5e-06, 11.2]

#Defino la cantidad de tuberías totales en el sistema
sist={'tub1':tub1,'tub2':tub2}

#gravedad m/s^2
g=9.81

# Asumir un valor inicial de f
f = 0.01

# Calcular f utilizando la iteración
error = 1e-6  # Tolerancia de error

Q= -1*(tub1[0]+tub2[0])#Caudal total

i=0
iterahf=[] #Inicializo las listas donde voy a almacenar los datos de las iteraciones
iterahe=[]
iteraQ=[]
iteraeffi=[]

for tub in sist.values():
    
    iteraeffi.append(tub[3]) #Almacena las eficiencias de cada tubería

    A = ((math.pi)/4)*(tub[4]**2) #Area de cada tuberia

    if tub == sist['tub1']:
        Q=Q
    
    else:
        l=str(i)
        x='tub'+l
        Q=Q+sist[x][0]

    V=Q/A

    Re = (V * tub[4]) / nu

    while True: #Entra para hallar el factor de friccion de la tubería

        f_new = (-2 * math.log10((tub[6]/ (3.7 * tub[4])) + (2.52 / (Re * math.sqrt(f))))) ** (-2)

        #V_new = math.sqrt((2*g*((E2-E1)+tub[2]))/((f_new*(tub[4]/tub[3])+tub[6])))
        #Re_new = (V * tub[4]) / nu
            
        print(f"el factor de fric (fc)) calculado es: {f_new:.6f}")

        # Verificar la convergencia con el error
        if abs(f - f_new) <= error:
            break 

        f = f_new
        #V = V_new
        #Re = Re_new

    hf=f*(tub[5]/tub[4])*((V**2)/(2*g)) #Calcula las perdidas por fricción de la tub
    he= tub[7]*((V**2)/(2*g)) #Calcula las perdidas menores de la tub

    print(f"el valor Q es: {Q:.6f}")
    print(f"el valor hf es: {V:.6f}")
    print(f"el valor hf es: {hf:.6f}")
    print(f"el valor he es: {he:.6f}")


    i=i+1

    iterahf.append(hf) #Alamacena los datos de las perdidas por fricción de cada iteración
    iterahe.append(he) #Alamacena los datos de las perdidas menores de cada iteración
    iteraQ.append(Q) #Almacena los datos de los caudales de cada iteración

hb=sum(iterahf)+sum(iterahe)-(E1-E2) #Calcula el el valor de la cabeza de la bomba teniendo en cuenta todos los valores calculados de los factores
print(f"el valor hb es: {hb:.6f} (m)")

for buscador, valor in enumerate(iteraeffi): #Variable que recorrerá la lista
    if valor !=0: #Condición donde la eficienca no es 0
        posicion = buscador
        efi=iteraeffi[posicion] #Asigna al valor de eficiencia el valor encontrado en la lista
    else:
        posicion=None
    break

if posicion is not None: #Si hay un valor de eficiencia va a calcular la potencia      
    Q=iteraQ[posicion] #El valor de caudal se extraerá de la misma posición donde encontró la eficienca, asi se asegura que la bomba esté en la tubería correcta.
    P=(rho*g)*Q*hb/efi #Se calcula finalmente la potencia con ese valor de Q y de efi.
    print(f"La potencia de la bomba será: {P:.6f} (W)")
else:
    print("No se encontraron valores de eficienca para calcular la bomba")
    


