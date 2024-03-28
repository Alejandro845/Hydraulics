import math

######ALGORITMO PARA EL CALCULO DEL FACTOR DE FRICCION CON PUNTO FIJO#####

# Valores de la tubería

epsilon = 0.00005  # Rugosidad de la tubería en metros
rho = 1000  # Densidad del fluido en kg/m^3
nu = 1.0e-6  # Viscosidad cinemática en m^2/s
V = 2.0  # Velocidad en m/s
D = 0.1  # Diámetro de la tubería en metros

# Calcular el número de Reynolds (Re)
Re = (V * D) / nu

# Verificar si el flujo es laminar o turbulento
if Re < 2000:
    # Flujo laminar
    f = Re / 64
else:
    # Asumir un valor inicial de f
    f = 0.01

    # Calcular f utilizando la iteración
    error = 1e-6  # Tolerancia de error

    while True:
        f_new = (-2 * math.log((epsilon / (3.7 * D)) + (2.52 / (Re * math.sqrt(f))))) ** -2

        # Verificar la convergencia
        if abs(f - f_new) <= error:
            break
        f = f_new

# Imprimir el factor de fricción calculado
print(f"El factor de fricción (f) calculado es: {f:.6f}")