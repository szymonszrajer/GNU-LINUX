import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt

N = 37848924 #populacja
I = 24687 # osoby zarażone na w pierwszym dniu badania
beta = 1.4 # liczba kontaktów dziennie
K = 14 # średni czas trwania infekcji wyrażony w dniach
time = 8 # liczba dni oznacza czas brany pod uwagę podzczas badań

R = 0 # początkowa liczba osób, które wyzdrowiały
S = N - I - R # liczba osób podatnych na zakażenie
k = 1/K # współczynnik wyzdrowień
t = np.linspace(0,time,time) # siatka punktów czasowych, oś x

# funkcja generująca równania różniczkowe
def derivatives(y, t, N, beta, k):
    S, I, R = y
    dSdt = -beta * S * I / N
    dIdt = beta * S * I / N - k * I
    dRdt = k * I
    return dSdt, dIdt, dRdt

# warunki początkowe wektora
y0 = S, I, R
# rozwiązanie równań różniczkowych za pomocą funkcji odeint
result = odeint(derivatives, y0, t, args=(N, beta, k))
# przypisanie wartości
S, I, R = result.T

# wykreślenie danych na trzech osobnych krzywych dla S(t), I(t) i R(t)
# kolor tła
backg = plt.figure(facecolor='w')
# opis osi
axis = backg.add_subplot(111, facecolor='#E5E8f9', axisbelow=True)
axis.set_xlabel('Time [days]')
axis.set_ylabel('Population (in 1 000 000s)')
# wykresy SIR, czas, wartości, kolor, szerokość, nazwa
axis.plot(t, S/1000000, '#46F0FF', lw=2, label='Susceptible')
axis.plot(t, I/1000000, '#FF7996', lw=2, label='Infected')
axis.plot(t, R/1000000, '#8EFF3B', lw=2, label='Recovered with immunity')
# ograniczenie osi y
axis.set_ylim(0, 38)
# "linie" - widoczne, które, kolor, szerokość, styl
axis.grid(b=True, which='both', c='w', lw=1, ls='-')
#legenda
legend = axis.legend()
plt.show()
print(I)