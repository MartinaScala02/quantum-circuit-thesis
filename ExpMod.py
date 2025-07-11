from qiskit import QuantumCircuit, QuantumRegister, ClassicalRegister, transpile
from qiskit_aer import AerSimulator
from qiskit.visualization import plot_histogram
import matplotlib.pyplot as plt


"-------------------------FUNZIONI-----------------------------"

def Sum():

    a = QuantumRegister(1)
    b = QuantumRegister(1)
    ci = QuantumRegister(1)

    qc = QuantumCircuit(ci, a, b, name = "Sum")
    

    qc.cx(ci[0], b[0])
    qc.cx(a[0], b[0])

    return qc.to_gate()



def Carry(): #opera sui singoli bit !!!

    a = QuantumRegister(1)
    b = QuantumRegister(1)
    c = QuantumRegister(1)
    c1 = QuantumRegister(1)
    
    qc = QuantumCircuit(c, a, b, c1, name = "Carry")


    qc.ccx(a[0], b[0], c1[0])
    qc.cx(a[0], b[0])
    qc.ccx(c[0], b[0], c1[0])
    qc.cx(a[0], b[0])

    return qc.to_gate()




def Add(n):


    a = QuantumRegister(n)
    b = QuantumRegister(n+1)
    c = QuantumRegister(n)
    qc = QuantumCircuit(c, a, b, name="Add")

    sum_gate = Sum()
    carry_gate = Carry()
    carry_gate_i = carry_gate.inverse(n)

    if n == 1:
        qc.append(carry_gate, [c[0], a[0], b[0], b[1]])
        qc.append(sum_gate, [c[0], a[0], b[0]])
    else:

        qc.append(carry_gate, [c[0], a[0], b[0], c[1]])

        add_ricorsivo = Add(n-1)
        qc.append(add_ricorsivo, c[1:n] + a[1:n] + b[1:n+1]) 

        qc.append(carry_gate_i, [c[0], a[0], b[0], c[1]])
        qc.append(sum_gate, [c[0], a[0], b[0]])

    return qc.decompose().to_gate()




  
def Copy(n): #copia il contenuto di un registro a in un altro registro b (b DEVE ESSERE INIZIALIZZATO A ZERO)
 
    a = QuantumRegister(n)
    b = QuantumRegister(n)
    qc = QuantumCircuit(a, b, name="Copy")


    for i in range(0, n):
        qc.cx(a[i], b[i])
        
    return qc.decompose().to_gate()




def Shift(n):

    a = QuantumRegister(n)
    qc = QuantumCircuit(a, name="Shift")
    
    for i in reversed(range(n-1)): #shifta i bit di 1 ciclicamente tramite swap
        qc.swap(a[i], a[i+1])
    
    return qc.decompose().to_gate()



def AddMod(n):


    a = QuantumRegister(n)
    b = QuantumRegister(n+1)
    c = QuantumRegister(n)
    M = QuantumRegister(n)

    t = QuantumRegister(1) #registro temporaneo
    

    qc = QuantumCircuit(c, a, b, M, t, name="AddMod")

    add_gate = Add(n)
    add_gate_c = Add(n).control(1)
    add_gate_i = add_gate.inverse(n)
    
    qc.append(add_gate, c[:] + a[:] + b[:])
    qc.append(add_gate_i, c[:] + M[:] + b[:])

    qc.cx(b[n], t[0])  

    qc.append(add_gate_c, [t[0]] + c[:] + M[:] + b[:])
    qc.append(add_gate_i, c[:] + a[:] + b[:])
    
    qc.x(b[n]) #voglio che il control not si attivi solo se b[n] è zero
    qc.cx(b[n], t[0])
    qc.x(b[n]) #lo faccio tornare normale
    
    qc.append(add_gate, c[:] + a[:] + b[:])


    return qc.decompose().to_gate()





def TimesMod(n, k):

    t = QuantumRegister(k) 
    a = QuantumRegister(n+1)
    b = QuantumRegister(k)
    c = QuantumRegister(n)
    M = QuantumRegister(n)
    p = QuantumRegister(n+1)


    t_a = QuantumRegister(n)
    t_b = QuantumRegister(1) #registro temporaneo per addmod

    
    qc = QuantumCircuit(c, a, b, M, p, t, t_a, t_b, name="TimesMod")

    add_gate = Add(n)
    add_gate_c = add_gate.control(1)
    add_gate_i = add_gate.inverse()
    addmod_gate = AddMod(n)
    addmod_gate_c = addmod_gate.control(1)
    add_gate_i_c = add_gate_i.control(1)   
    shift_gate = Shift(n+1)
    shift_gate_i = shift_gate.inverse()

    

    for i in range(0, k):
        qc.append(add_gate_i, c[:] + M[:] + a[:])
        qc.cx(a[n], t[i])
        qc.append(add_gate_c, [t[i]] + c[:] + M[:] + a[:])
        qc.append(addmod_gate_c, [b[i]] + t_a[:] + a[0:n] + p[:] + M[:] + t_b[:]) 
        qc.append(shift_gate, a[:])
 

    for i in reversed(range(0, k)):
        qc.append(shift_gate_i, a[:])
        qc.append(add_gate_i_c, [t[i]] + c[:] + M[:] + a[:])
        qc.cx(a[n], t[i])
        qc.append(add_gate, c[:] + M[:] + a[:])

    return qc.decompose().to_gate()



def SquareMod(n):

    a = QuantumRegister(n+1)
    M = QuantumRegister(n)
    s = QuantumRegister(n+1)

    c = QuantumRegister(n) 
    t = QuantumRegister(n)
    t1 = QuantumRegister(n)
    t_a = QuantumRegister(n) 
    t_b = QuantumRegister(1) 


    qc = QuantumCircuit(c, a, M, s, t, t1, t_a, t_b, name= "SquareMod")

    copy_gate = Copy(n)
    copy_gate_i = copy_gate.inverse()
    timesmod_gate = TimesMod(n, n)
    
    qc.append(copy_gate, a[0:n] + t[:])
    qc.append(timesmod_gate, c[:] + a[:] + t[:] + M[:] + s[:] + t1[:] + t_a[:] + t_b[:])
    qc.append(copy_gate_i, a[0:n] + t[:])

    return qc.decompose().to_gate()



def ExpMod(n, k):


    a = QuantumRegister(n+1)
    b = QuantumRegister(k)
    M = QuantumRegister(n)
    p = QuantumRegister(n+1)
    e = QuantumRegister(n+1)
   
    c = QuantumRegister(n) 
    t = QuantumRegister(n+1)
    t_a = QuantumRegister(n) 
    t_b = QuantumRegister(1)


    copy_gate = Copy(n+1)
    copy_gate_c = copy_gate.control(1)
    copy_gate_i = copy_gate.inverse()
    copy_gate_i_c = copy_gate_i.control(1)
    timesmod_gate = TimesMod(n, n+1)
    timesmod_gate_i = timesmod_gate.inverse()
    timesmod_gate_c = timesmod_gate.control(1)
    timesmod_gate_i_c = timesmod_gate_i.control(1)
    squaremod_gate = SquareMod(n)
    squaremod_gate_i = squaremod_gate.inverse()


    if k == 1:
        qc = QuantumCircuit(c, a, b, M, p, e, t, t_a, t_b, name= "ExpMod")
        qc.x(b[0]) #lo nego
        qc.append(copy_gate_c, [b[0]] + p[:] + e[:])
        qc.x(b[0]) #lo faccio tornare normale
        qc.append(timesmod_gate_c, [b[0]] + c[:] + a[:] + p[:] + M[:] + e[:] + t[:] + t_a[:] + t_b[:]) 

    else:

        u = QuantumRegister((k-1)*(n+1)) #registro temporaneo expmod
        v = QuantumRegister((k-1)*(n+1)) #registro temporaneo expmod

        ind1 = (k-2)*(n+1)
        ind2 = (k-1)*(n+1)

        qc = QuantumCircuit(c, a, b, M, p, e, t, t_a, t_b, u, v, name= "ExpMod")

        qc.x(b[0])
        qc.append(copy_gate_c, [b[0]] + p[:] + v[:])
        qc.x(b[0]) 
        qc.append(timesmod_gate_c, [b[0]] + c[:] + a[:] + p[:] + M[:] + e[:] + t[:] + t_a[:] + t_b[:]) 
        qc.append(squaremod_gate, c[:] + a[:] + M[:] + u[ind1:ind2] + v[0:n] + t[0:n] + t_a[:] + t_b[:])
        expmod_gate_recursive = ExpMod(n, k-1)
        qc.append(expmod_gate_recursive, c[:] + u[ind1:ind2] + b[1:k] + M[:] + v[ind1:ind2] + e[:] + t[:] + t_a[:] + t_b[:] + u[0:ind1] + v[0:ind1]) 
        qc.append(squaremod_gate_i, c[:] + a[:] + M[:] + u[ind1:ind2] + v[0:n] + t[0:n] + t_a[:] + t_b[:])
        qc.append(timesmod_gate_i_c, [b[0]] + c[:] + a[:] + p[:] + M[:] + e[:] + t[:] + t_a[:] + t_b[:])

        qc.x(b[0])
        qc.append(copy_gate_i_c, [b[0]] + p[:] + v[:])
        qc.x(b[0]) 

    return qc.decompose().to_gate()

"-------------------------------------------------------------------"

n = 2
k = n

num_q=2*(k-1)*(n+1)+7*n+5+k
a = QuantumRegister(n+1)
b = QuantumRegister(k)
c = QuantumRegister(n)
p = QuantumRegister(n+1)
e = QuantumRegister(n+1)
u = QuantumRegister((k-1)*(n+1))
v = QuantumRegister((k-1)*(n+1))
M = QuantumRegister(n)
t = QuantumRegister(n+1)
t1 = QuantumRegister(n)
t_a = QuantumRegister(n)
t_b = QuantumRegister(1)

creg = ClassicalRegister(n+1)


qc = QuantumCircuit(a, b, c, M, p, e, t, t1, t_a, t_b, u, v, creg)



"-------------------INIZIALIZZAZIONE------------------------"


qc.x(a[0])
qc.x(a[1])
#qc.x(a[2])
#qc.x(a[3])
qc.x(b[0])
#qc.x(b[1])
#qc.x(b[2])
qc.x(p[0])
qc.x(M[0])
qc.x(M[1])
#qc.x(M[2])
#qc.x(M[3])

"-------------------------------------------------------------"

expmod_gate = ExpMod(n, k)
qc.append(expmod_gate, list(range(0,num_q)))

qc.measure(e, creg)
qc.draw(output="mpl")

simulator = AerSimulator()
qc1 = transpile(qc, backend=simulator)  

job = simulator.run(qc1, shots=1024)
counts = job.result().get_counts()
plot_histogram(counts)
plt.show()



