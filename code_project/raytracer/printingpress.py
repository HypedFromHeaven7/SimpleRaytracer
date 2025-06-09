import pyperclip as ppr

a=    [[0.1, 0.1, 0],
    [0, 0, 0],
    [-0.1, -0.1, 0]]

opa = "rays = ["
for n in range(len(a)):
    if n == 0:
        opa += f"Ray({str(a[n])}),"
    else:
        opa += f"\n            Ray({str(a[n])}),"

opa += "]"

print(opa)

ppr.copy(opa)
