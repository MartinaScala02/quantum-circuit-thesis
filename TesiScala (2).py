from qiskit import QuantumCircuit, QuantumRegister, ClassicalRegister, transpile
from qiskit_aer import AerSimulator
from qiskit.visualization import plot_histogram
import matplotlib.pyplot as plt

# Definizione porte
# Operate on three single qubits ordered as c,a,b
def Sum():
    a=QuantumRegister(1)
    b=QuantumRegister(1)
    c=QuantumRegister(1)
    A=QuantumCircuit(c,a,b,name="Sum")
    A.cx(a,b)
    A.cx(c,b)
    return A.decompose().to_gate()

# the carry C(a,b,c) is 1 if two or more of the bits a, b, c, are 1
# Operate on four single qubits
def Carry():
    a=QuantumRegister(1)
    b=QuantumRegister(1)
    c=QuantumRegister(1)
    c1=QuantumRegister(1)
    A=QuantumCircuit(c,a,b,c1,name="Carry")
    A.ccx(a,b,c1)
    A.cx(a,b)
    A.ccx(c,b,c1)
    A.cx(a,b)
    return A.decompose().to_gate()

#  adds two n-bit binary numbers
#  where a and c are n-qubit registers and b is an (n+1)-qubit register, adds two n-bit numbers,
#  placed in registers a and b, and puts the result in register b when register c and the highest order
#  bit, bn, of register b are initially 0.
def Add(n):
    a=QuantumRegister(n)
    b=QuantumRegister(n+1)
    c=QuantumRegister(n)

    A=QuantumCircuit(c,a,b, name="Add")

    CG=Carry()
    CG_i=CG.inverse(n)
    S=Sum()

    if(n==1):
        A.append(CG,c[:]+a[:]+[b[0]]+[b[1]])
        A.append(S,c[:]+a[:]+[b[0]])
    else:
        A.append(CG,[c[0]]+[a[0]]+[b[0]]+[c[1]])
        Ad1=Add(n-1)
        A.append(Ad1,c[1:n]+a[1:n]+b[1:n+1])
        A.append(CG_i,[c[0]]+[a[0]]+[b[0]]+[c[1]])
        A.append(S,[c[0]]+[a[0]]+[b[0]])
    return A.decompose().to_gate()


#  The following program defines modular addition for n-bit binary numbers a and b
# where the registers a and M have n qubits and b is an n+1-qubit register. When the highest order
# bit, bn, of register b is initially 0, the transformation AddMod replaces the contents of register
# b with b+amodM, where M is the contents of register M. The contents of registers a and M
# (and the temporaries c and t) are unchanged by AddMod.
def AddMod(n):
    t=QuantumRegister(1)
    c=QuantumRegister(n)
    a=QuantumRegister(n)
    b=QuantumRegister(n+1)
    M=QuantumRegister(n)

    A = QuantumCircuit(a,b,M,c,t)

    AG=Add(n)
    AG_i=AG.inverse(n)
    AG_c=AG.control(1)

    A.append(AG,c[:]+a[:]+b[:])
    A.append(AG_i,c[:]+M[:]+b[:])
    A.cx(b[n],t)
    A.append(AG_c,[t[0]]+c[:]+M[:]+b[:])
    A.append(AG_i,c[:]+a[:]+b[:])
    A.x(b[n])
    A.cx(b[n],t)
    A.x(b[n])
    A.append(AG,c[:]+a[:]+b[:])

    return A.decompose().to_gate()

def Shift(n):
    a = QuantumRegister(n)
    qc = QuantumCircuit(a, name="Shift")
    
    for i in reversed(range(n-1)): 
        qc.swap(a[i], a[i+1])
    
    return qc.decompose().to_gate()

def TimesMod(n,k):
    t=QuantumRegister(k)
    c=QuantumRegister(n)
    a=QuantumRegister(n+1)
    b=QuantumRegister(k)
    M=QuantumRegister(n)
    p=QuantumRegister(n+1)

    #Temporary register for AddMod
    t_a=QuantumRegister(n)
    t_b=QuantumRegister(1)

    A=QuantumCircuit(a,b,M,p,c,t,t_a,t_b)

    AG=Add(n)
    AMG=AddMod(n)
    AG_i=AG.inverse(n)
    AG_c=AG.control(1)
    AG_i_c=AG_i.control(1)
    AMG_c=AMG.control(1)
    SH=Shift(n+1)
    SH_i=SH.inverse(n+1)
    

    for i in range(0,k):
        A.append(AG_i,c[:]+M[:]+a[:])
        A.cx(a[n],t[i])
        A.append(AG_c,[t[i]]+c[:]+M[:]+a[:])
        A.append(AMG_c,[b[i]]+a[0:n]+p[:]+M[:]+t_a[:]+t_b[:])
        A.append(SH,a)

    for i in reversed(range(0,k)):
        A.append(SH_i,a)
        A.append(AG_i_c,[t[i]]+c[:]+M[:]+a[:])
        A.cx(a[n],t[i])
        A.append(AG,c[:]+M[:]+a[:])

    return A.decompose().to_gate()

def Copy(n):
    a=QuantumRegister(n)
    b=QuantumRegister(n)
    A=QuantumCircuit(a,b)
    
    for i in range(0,n):
        A.cx(a[i],b[i])

    return A.decompose().to_gate()

def SquareMod(n):
    t=QuantumRegister(n)
    a=QuantumRegister(n+1)
    M=QuantumRegister(n)
    s=QuantumRegister(n+1)
    

    c=QuantumRegister(n)
    t1=QuantumRegister(n)
    t_a=QuantumRegister(n)
    t_b=QuantumRegister(1)

    A=QuantumCircuit(a,M,s,t,c,t1,t_a,t_b)

    CG=Copy(n)
    CG_i=CG.inverse(n)
    TM=TimesMod(n,n)

    A.append(CG,a[0:n]+t[:])
    A.append(TM,a[:]+t[:]+M[:]+s[:]+c[:]+t1[:]+t_a[:]+t_b[:])
    A.append(CG_i,a[0:n]+t[:])

    return A.decompose().to_gate()


def ExpMod(n,k):
    # Registri computazionali
    a=QuantumRegister(n+1)
    b=QuantumRegister(k)
    M=QuantumRegister(n)
    p=QuantumRegister(n+1)
    e=QuantumRegister(n+1)

    # Registri temporanei
    c=QuantumRegister(n)
    t=QuantumRegister(n+1)
    t_a=QuantumRegister(n)
    t_b=QuantumRegister(1)


    CG=Copy(n+1)
    CG_i=CG.inverse()
    CG_c=CG.control(1)
    CG_i_c=CG_i.control(1)
    TMG=TimesMod(n,n+1)
    TMG_c=TMG.control(1)
    SMG=SquareMod(n)
    SMG_i=SMG.inverse()
    TMG_i=TMG.inverse()
    TMG_i_c=TMG_i.control(1)

    if(k==1):
        A=QuantumCircuit(a,b,M,p,e,c,t,t_a,t_b)
        A.x(b[0])
        A.append(CG_c,[b[0]]+p[:]+e[:])
        A.x(b[0])
        A.append(TMG_c,[b[0]]+a[:]+p[:]+M[:]+e[:]+c[:]+t[:]+t_a[:]+t_b[:])
    else:
        u=QuantumRegister((k-1)*(n+1))
        v=QuantumRegister((k-1)*(n+1))
        ind1=(k-2)*(n+1)
        ind2=(k-1)*(n+1)
        A=QuantumCircuit(a,b,M,p,e,c,t,t_a,t_b,u,v)
        A.x(b[0])
        A.append(CG,p[:]+v[:])
        A.x(b[0])
        A.append(TMG_c,[b[0]]+a[:]+p[:]+M[:]+e[:]+c[:]+t[:]+t_a[:]+t_b[:])
        A.append(SMG,a[:]+M[:]+u[ind1:ind2]+v[0:n]+c[:]+t[0:n]+t_a[:]+t_b[:]) # v come registro temporaneo tanto viene resettato
        Exp1=ExpMod(n,k-1)
        A.append(Exp1,u[ind1:ind2]+b[1:k]+M[:]+v[ind1:ind2]+e[:]+c[:]+t[:]+t_a[:]+t_b[:]+u[0:ind1]+v[0:ind1])
        A.append(SMG_i,a[:]+M[:]+u[ind1:ind2]+v[0:n]+c[:]+t[0:n]+t_a[:]+t_b[:])
        A.append(TMG_i_c,[b[0]]+a[:]+p[:]+M[:]+e[:]+c[:]+t[:]+t_a[:]+t_b[:])
        A.x(b[0])
        A.append(CG_i_c,[b[0]]+p[:]+v[:])
        A.x(b[0])
    return A.decompose().to_gate()




n=2
k=2
num_q=2*(k-1)*(n+1)+7*n+5+k
a=QuantumRegister(n+1)
b=QuantumRegister(k)
M=QuantumRegister(n)
p=QuantumRegister(n+1)
e=QuantumRegister(n+1)
c=QuantumRegister(n)
t=QuantumRegister(n+1)
t1=QuantumRegister(n)
t_a=QuantumRegister(n)
t_b=QuantumRegister(1)
u=QuantumRegister((k-1)*(n+1))
v=QuantumRegister((k-1)*(n+1))
res1=ClassicalRegister(n+1)
qc=QuantumCircuit(a,b,M,p,e,c,t,t_a,t_b,u,v,res1)

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
SG=ExpMod(n,k)
print(num_q)
qc.append(SG,list(range(0,num_q)))
qc.measure(e,res1)
qc.draw(output="mpl")

simulator = AerSimulator()
qc1 = transpile(qc, backend=simulator)  

job = simulator.run(qc1, shots=1024)
counts = job.result().get_counts()
plot_histogram(counts)
plt.show()



