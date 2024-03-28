import math

########## CALCULO DE COMPROBACIÓN DE DISEÑO PARA TUBERÍÁS SIMPLES ##########

#Valores de entrada

epsilon = 1.5e-06 #Rugosidad de la tubería en metros
rho = 998.2 # Densidad del fluido en kg/m^3
mu = 0.001005  #Viscosidad dinámica 
nu = mu/rho # Viscosidad cinemática en m^2/s


###Bombas###    
hb = 0  #Cabeza de la bomba
ht = 0  #Cabeza de la turbina


D = 0.293  # Diámetro de la tubería en metros
E1 = 43.5  #Energia en 1
E2 = 0 #Energia en 2
L = 730 #longitud de la tuberia en metros
K = 11.8  #Sumatoria de coeficiente de perdidas menores

#gravedad m/s^2
g=9.81

# Asumir un valor inicial de f
f = 0.01

# Calcular f utilizando la iteración
error = 1e-6  # Tolerancia de error

V=math.sqrt((2*g*(E1-E2+hb+ht))*(((f*(L/D))+K)**(-1.)))

Re = (V * D) / nu


while True:
    	
    t1=epsilon/(3.7*D)
    t2=2.51/(Re*math.sqrt(f))

    #f_new = (-2 * math.log10((epsilon / (3.7 * D)) + (2.52 / (Re * math.sqrt(f))))) ** (-2)   

    f_new=(-2.*math.log10(t1+t2))**(-2.)

    V_new = math.sqrt((2*g*(E1-E2+hb+ht))/((f_new*(L/D)+K)))

    Re_new = (V_new * D) / nu

    print(f"el factor de fric (fc)) calculada es: {f_new:.6f}")
    print(f"La velocidad (Vc)) calculada es: {V:.6f}")

        # Verificar la convergencia
    if abs(f - f_new) <= error:
        break 

    f = f_new
    V = V_new
    Re = Re_new

# Imprimir el factor de fricción calculado
print(f"El factor de fricción (f) calculado es: {f:.6f}")
# Imprimir la velocidad calculada
print(f"La velocidad (V)) calculada es: {V:.6f}")









































