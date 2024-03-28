import math

########## CALCULO DE POTENCIA PARA TUBERÍÁS SIMPLES #####

#Valores de entrada#

epsilon = 4.6e-05 #Rugosidad de la tubería en metros
rho = 999.1 # Densidad del fluido en kg/m^3
mu = 0.001005  #Viscosidad dinámica 
#nu = mu/rho # Viscosidad cinemática en m^2/s (Si se necesita calcular)
nu = 1.14e-06 # Viscosidad cinemática en m^2/s
Q = 0.46 #Caudal
D = 0.3556  # Diámetro de la tubería en metros
E1 = 0  #Energia en 1 m
E2 = 15  #Energia en 2 m
L = 370 #longitud de la tuberia en metros
K = 7.4  #Sumatoria de coeficiente de perdidas menores

############
###Bombas###    
############

hb = 0  #Cabeza de la bomba m
ht = 0  #Cabeza de la turbina m
efi = 0.75  #Eficiencia de la bomba/turbina

 
g=9.81 #gravedad m/s^2

# Asumir un valor inicial de f
f = 0.01

# Calcular f utilizando la iteración
error = 1e-6  # Tolerancia de error

# Calculo del area y la velocidad
A = ((math.pi)/4)*(D**2)
V=Q/A

# Calculo del numero de Reynolds
Re = (V * D) / nu


while True:

    f_new = (-2 * math.log10((epsilon / (3.7 * D)) + (2.52 / (Re * math.sqrt(f))))) ** (-2)

    print(f"el factor de fric (fc)) calculada es: {f_new:.6f}")

    # Verificar la convergencia con el error
    if abs(f - f_new) <= error:
        break 

    f = f_new

 #Calculo de las perdidas por friccion con Darcy-Weischback
    #Calculo de las perdidas menores

    #Calculo de la potencia 

    #Calculo de la cabeza de la bom
hf= f*(L/D)*((V**2)/(2*g))
he= K*((V**2)/(2*g))
hb=hf+he+(E2-E1)
P=(rho*g)*Q*hb/efi


# Imprimir el factor de fricción calculado
print(f"El factor de fricción (f) calculado es: {f:.6f}")

# Imprimir las perdidas por friccion
print(f"El factor de fricción (hf) calculado es: {hf:.6f}")

# Imprimir las perdidas menores
print(f"El factor de fricción (he) calculado es: {he:.6f}")

# Imprimir la cabeza de bomba
print(f"El factor de fricción (hb) calculado es: {hb:.6f}")

# Imprimir la potencia
print(f"El factor de fricción (P) calculado es: {P:.6f}")


