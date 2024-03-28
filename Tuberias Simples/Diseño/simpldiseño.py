import math

########## CALCULO DE DISEÑO PARA TUBERÍAS SIMPLES ##########

#Valores de entrada

epsilon = 0.00015 #Rugosidad de la tubería en metros
rho = 999.3 # Densidad del fluido en kg/m^3
mu = 0.00117  #Viscosidad dinámica 
nu = mu/rho # Viscosidad cinemática en m^2/s
Q = 0.12 #Caudal teorico a spuerar o alcanzar


###Bombas###    
hb = 0  #Cabeza de la bomba
ht = 0  #Cabeza de la turbina
efi = 0  #Eficiencia de la bomba/turbina


E1 = 2.2  #Energia en 1
E2 = 0 #Energia en 2
L = 150 #longitud de la tuberia en metros
K = 0.5+0.8+1+1 #Sumatoria de coeficiente de perdidas menores

#gravedad m/s^2
g=9.81

# Asumir un valor inicial de f
f = 0.01

# Tolerancia de error
error = 1e-6  


#Diametros comerciales en pulgadas/mm

CD={'3/4':23.63,'1':30.20,'1 1/4':38.14,'1 1/2':43.68,'2':54.58,'2 1/12':66.07,'3':80.42,'4':103.42,'6':152.22,'8':198.21,'10':247.09,'12':293.07,'14':321.76,'16':367.70,'18':413.66,'20':459.64,'24':551.54}

#Creo una lista para almacenar solo los valores de los diametros
diametros=[]

#Ciclo para extraer y almcanear los valores en la lista de diametros
for diam in CD:
    diametros.append(float(CD[diam]))


#Ciclo que recorre cada uno de los diametros comerciales de la lista
for i in range(len(diametros)):
    D=diametros[i]

    D=D/1000

    V=math.sqrt((2*g*(E1-E2+hb+ht))/((f*(L/D))+K))

    Re = (V * D) / nu

    #While para calcular el valor de punto fijo con punto fijo
    while True:

        f_new = (-2 * math.log10((epsilon / (3.7 * D)) + (2.52 / (Re * math.sqrt(f))))) ** (-2)   

        V_new = math.sqrt((2*g*(E1-E2+hb+ht))/((f_new*(L/D)+K)))

        Re_new = (V * D) / nu
 
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

Q=Qnew
hf= f*(L/D)*((V**2)/(2*g))
he= K*((V**2)/(2*g))


# Imprimir el factor de fricción calculado
print(f"El factor de fricción (f) calculado es: {f:.6f}")
# Imprimir la velocidad calculada
print(f"La velocidad (V)) calculada es: {V:.6f}")
# Imprimir el diametro con el que sale
print(f"El diametro (D)) calculado es: {D:.6f}")
# Imprimir el caudal con el que sale
print(f"El caudal (Q)) calculado es: {Q:.6f}")
# Imprimir las perdidas por friccion
print(f"El factor de fricción (hf) calculado es: {hf:.6f}")
# Imprimir las perdidas menores
print(f"Las perdidas menores (he) calculadas son: {he:.6f}")




