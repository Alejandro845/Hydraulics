import math

########## CALCULO DE COMPROBACIÓN DE DISEÑO PARA TUBERÍAS EN SERIE - EJEMPLO 2 ##########

#Datos globales del sistema

rho = 999.3 # Densidad del fluido en kg/m^3
mu = 0.001005  #Viscosidad dinámica 
#nu = mu/rho # Viscosidad cinemática en m^2/s
nu = 1e-06 # Viscosidad cinemática en m^2/s
ki = 0 #Coeficiente de perdidas menores en la entrada
ko = 0 #Coeficiente de perdidas menores en la salida
E1 = 10  #Energia en 1
E2 = 0 #Energia en 2


# Datos locales de cada una de las tuberias

#tub_n=[Qo (Caudal en la salida o lateral), Qi (Caudal de entrada de la tuberia), hb o ht (cabeza de bomba o turb), D (diametro), L (long), epsilon (rugosidad), K (sumatoria de perdidas menores)]
#tub_n=[Qo, Qi, hb/ht, D, L, epsilon, K]

tub1=[0, 0, 0, 0.3, 300, 0.00025, 0]
tub2=[0, 0, 0, 0.5, 300, 0.00025, 0]
tub3=[0, 0, 0, 0.3, 600, 0.00025, 0]
#tub4=[0, 0, 0, 0.254, 121, 1.5e-06, 7.5]
sist={'tub1':tub1,'tub2':tub2,'tub3':tub3} #'tub4':tub4 


#gravedad m/s^2
g=9.81

# Asumir un valor inicial de f
f = 0.01

# Calcular f utilizando la iteración
error = 1e-6  # Tolerancia de error

HT=E1-E2

j=0
delta=0
while True:
    i=0
    for tub in sist.values():
        A = ((math.pi)/4)*(tub[3]**2)
        

        
        if tub == sist['tub1']: #Solo corrige el valor de hf si esta en la primera tuberia
            denohf=0
            a=0
            for tuba in sist.values(): #tuba solo itera para obtener la sumatoria de cada una de las longitudes y diametros de las tuberias
                a=(tuba[4])/((tuba[3])**5)
                denohf+=a #sumatoria

            hf=HT*(((tub[4])/((tub[3])**5))/(denohf))
           
            if j>=1 and i==0:
                aaa=(HT-HG)*(((tub[4])/((tub[3])**5))/(denohf))
                delta+=aaa
                hf=hf+delta 
            HG=0
            V=((-2*math.sqrt(2*g*tub[3]*hf))/(math.sqrt(tub[4])))*(math.log10(((tub[5])/(3.7*tub[3]))+((2.52*nu*math.sqrt(tub[4]))/(tub[3]*math.sqrt(2*g*tub[3]*hf)))))
            
        
            Q=V*A  
        else:  #string para restar los caudales laterales de cada tuberia
            l=str(i)
            x='tub'+l
            Q=Q+sist[x][0]
            V=Q/A
        
        he= tub[6]*((V**2)/(2*g))

        Re = (V * tub[3]) / nu

        while True:

            f_new = (-2 * math.log10((tub[5]/ (3.7 * tub[3])) + (2.52 / (Re * math.sqrt(f))))) ** (-2)

            #V_new = math.sqrt((2*g*((hf+he)+tub[2]))/((f_new*(tub[4]/tub[3])+tub[6])))
            #Re_new = (V * tub[3]) / nu
            
            print(f"el factor de fric (fc)) calculada es: {f_new:.6f}")

            # Verificar la convergencia con el error
            if abs(f - f_new) <= error:
                break 

            f = f_new
            #V = V_new
            #Re = Re_new

        print(f"el caudal (Q) al final de la iteracion es: {Q:.6f}")

        he= tub[6]*((V**2)/(2*g))
        hf=f*(tub[4]/tub[3])*((V**2)/(2*g))
        if tub == sist['tub3']: #valor maximo de tuberias
            HG+=ko*((V**2)/(2*g))+he+hf
        if tub == sist['tub1']:
            HG+=ki*((V**2)/(2*g))+he+hf
        if tub == sist['tub2']:
            HG+=he+hf

        print(f"el valor hf es: {hf:.6f}")
        print(f"el valor he es: {he:.6f}")
        print(f"el valor de HG es: {HG:.6f}")
        i=i+1
    
    if abs(HG-HT)<=error:
        break

    j=j+1
    

