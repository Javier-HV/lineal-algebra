
from sympy import symbols, simplify, solve, Mul, Matrix, S, sqrt, zeros, eye, diag, latex, parse_expr, Rational, Eq
from sympy import *
import numpy as np
from IPython.display import display, Markdown


# imprimir_operacion: Función para imprimir una operación y su resultado en formato Markdown y guardarlo en una lista.
def imprimir_operacion(operacion, A, output_content):
# Parámetros: operacion (str), A (Matrix) - Imprime la operación y la matriz A

    display(Markdown(f"Operación: ${operacion}$"))
    display(Matrix(A))
    output_content.append(f"Operación: ${operacion}$\n")
    output_content.append(f"{latex(Matrix(A))}\n")
    
    return output_content

    
# rrefsymb: Función para calcular la forma escalonada reducida por filas de una matriz, imprimirla y guardarla en un archivo de texto.
def rrefsymb(M, archivo_salida="RREF.txt"):
# Parámetros: M (Matrix), archivo_salida (str) - Calcula RREF de M y lo guarda en archivo_salida

    A = Matrix(M)
    n, m = A.shape
    output_content = []  # Lista para almacenar el contenido impreso

    fila_pivote = 0
    for i in range(m):
        # Encontrar un pivote no cero en la columna actual o en las siguientes
        pivote_encontrado = False
        for j in range(fila_pivote, n):
            if A[j, i] != 0:
                pivote_encontrado = True
                # Mover la fila con el pivote a la posición fila_pivote
                if j != fila_pivote:
                    A.row_swap(j, fila_pivote)
                    output_content = imprimir_operacion(f"R_{{{j + 1}}} \\leftrightarrow R_{{{fila_pivote + 1}}}", A, output_content)
                break
        if not pivote_encontrado:
            continue  # No se encontró un pivote en esta columna, continuar a la siguiente

        # Hacer el pivote igual a 1
        if A[fila_pivote, i] != 1:
            factor = S(1) / A[fila_pivote, i]
            A.row_op(fila_pivote, lambda v, j: v * factor)
            output_content= imprimir_operacion(f"\\frac{{1}}{{{A[fila_pivote, i]}}} R_{{{fila_pivote + 1}}} \\rightarrow R_{{{fila_pivote + 1}}}", A, output_content)

        # Hacer ceros en todas las demás filas de esta columna
        for j in range(n):
            if j != fila_pivote and A[j, i] != 0:
                factor = -A[j, i]
                A.row_op(j, lambda v, k: v + factor * A[fila_pivote, k])
                output_content = imprimir_operacion(f"R_{{{j + 1}}} + ({-factor}) R_{{{fila_pivote + 1}}} \\rightarrow R_{{{j + 1}}}", A, output_content)

        fila_pivote += 1

    display(Markdown("Forma escalonada reducida (RREF):"))
    display(Matrix(A))
    output_content.append("Forma escalonada reducida (RREF):\n")
    output_content.append(f"{latex(Matrix(A))}\n")

    with open(archivo_salida, "w", encoding='utf-8') as file:
        file.writelines(output_content)

    return A

# reduccion_triangular: Descripción no disponible.
def reduccion_triangular(M, archivo_salida="reduccion_triangular.txt"):



    A = Matrix(M)
    n, m = A.shape
    output_content = []  # Lista para almacenar el contenido impreso

    fila_pivote = 0
    for i in range(m):
        # Encontrar un pivote no cero en la columna actual o en las siguientes
        pivote_encontrado = False
        for j in range(fila_pivote, n):
            if A[j, i] != 0:
                pivote_encontrado = True
                # Mover la fila con el pivote a la posición fila_pivote
                if j != fila_pivote:
                    A.row_swap(j, fila_pivote)
                    output_content= imprimir_operacion(f"R_{{{j + 1}}} \\leftrightarrow R_{{{fila_pivote + 1}}}", A, output_content)
                break

        if not pivote_encontrado:
            continue  # No se encontró un pivote en esta columna, continuar a la siguiente

        # Hacer el pivote igual a 1
        if A[fila_pivote, i] != 1:
            factor = S(1) / A[fila_pivote, i]
            A.row_op(fila_pivote, lambda v, j: v * factor)
            output_content = imprimir_operacion(f"\\frac{{1}}{{{A[fila_pivote, i]}}} R_{{{fila_pivote + 1}}} \\rightarrow R_{{{fila_pivote + 1}}}", A, output_content)

        # Hacer ceros debajo del pivote
        for j in range(fila_pivote + 1, n):
            if A[j, i] != 0:
                factor = -A[j, i]
                A.row_op(j, lambda v, k: v + factor * A[fila_pivote, k])
                output_content = imprimir_operacion(f"R_{{{j + 1}}} + ({-factor}) R_{{{fila_pivote + 1}}} \\rightarrow R_{{{j + 1}}}", A, output_content)

        fila_pivote += 1

    display(Markdown("Matriz triangular superior:"))
    display(Matrix(A))
    output_content.append("Matriz triangular superior:\n")
    output_content.append(f"{latex(Matrix(A))}\n")

    with open(archivo_salida, "w", encoding='utf-8') as file:
        file.writelines(output_content)

    return A

# inversa: Descripción no disponible.
def inversa(M, archivo_salida="INVERSA.txt"):
    A = Matrix(M)
    n, m = A.shape

    output_content = []  # Lista para almacenar todo el contenido que se imprime y guarda

    if n != m:
        display(Markdown("La matriz debe ser cuadrada"))
        output_content.append("La matriz debe ser cuadrada\n")
        return

    # Crear la matriz aumentada con la matriz identidad
    matriz_identidad = eye(n)
    matriz_aumentada = A.row_join(matriz_identidad)

    display(Markdown("Matriz aumentada inicial:"))
    display(matriz_aumentada)
    output_content.append("Matriz aumentada inicial:\n")
    output_content.append(f"{latex(matriz_aumentada)}\n")

    # Aplicar operaciones de fila
    for i in range(n):
        # Hacer el elemento diagonal [i, i] igual a 1
        if matriz_aumentada[i, i] == 0:
            display(Markdown(f"Intercambio de fila necesario en la Fila ${i+1}$"))
            output_content.append(f"Intercambio de fila necesario en la Fila ${i+1}$\n")
            for j in range(i+1, n):
                if matriz_aumentada[j, i] != 0:
                    matriz_aumentada.row_swap(i, j)
                    display(matriz_aumentada)
                    output_content.append(f"{latex(matriz_aumentada)}\n")
                    break

        # Hacer el elemento diagonal [i, i] igual a 1 si no lo es
        if matriz_aumentada[i, i] != 1:
            factor = 1 / matriz_aumentada[i, i]
            matriz_aumentada.row_op(i, lambda v, j: v * factor)
            display(Markdown(f"Haciendo el elemento diagonal [{i+1}, {i+1}] igual a 1"))
            display(matriz_aumentada)
            output_content.append(f"Haciendo el elemento diagonal [{i+1}, {i+1}] igual a 1\n")
            output_content.append(f"{latex(matriz_aumentada)}\n")

        # Hacer ceros arriba y debajo del elemento diagonal [i, i]
        for j in range(n):
            if j != i and matriz_aumentada[j, i] != 0:
                factor = -matriz_aumentada[j, i]
                matriz_aumentada.row_op(j, lambda v, k: v + factor * matriz_aumentada[i, k])
                display(Markdown(f"Haciendo ceros en la columna {i+1} excepto en la Fila ${i+1}$"))
                display(matriz_aumentada)
                output_content.append(f"Haciendo ceros en la columna {i+1} excepto en la Fila ${i+1}$\n")
                output_content.append(f"{latex(matriz_aumentada)}\n")

    # Separar la matriz inversa
    inversa = matriz_aumentada[:, n:]

    display(Markdown("Matriz inversa:"))
    display(inversa)
    output_content.append("Matriz inversa:\n")
    output_content.append(f"{latex(inversa)}\n")

    # Guardar todo el contenido en un archivo en la carpeta raíz
    with open(archivo_salida, "w", encoding='utf-8') as file:
        file.writelines(output_content)

    return inversa

# solve_GJ: Descripción no disponible.
def solve_GJ(A, b, archivo_salida="solucion_sistema.txt", triang=False):

    # Crear la matriz ampliada con A y b
    matriz_ampliada = Matrix(A).row_join(Matrix(b))
    n, m = matriz_ampliada.shape

    output_content = []
    if triang:
        A_triangular = Matrix(A)
        b_triangular = Matrix(b)
    else:
        # Aplicar la reducción a matriz triangular
        display(Markdown("---"))
        matriz_triangular = reduccion_triangular(matriz_ampliada, archivo_salida="TRIANGULAR.txt")
        display(Markdown("---"))

        A_triangular = matriz_triangular[:, :-1]
        b_triangular = matriz_triangular[:, -1]

    Esol = True
    inf_sol = False
    columnas_cero_indices = [j for j in range(m-1) if all(A_triangular[:, j][i] == 0 for i in range(n))]

    for i in range(n):
        if all(A_triangular.row(i)[j] == 0 for j in range(m-1)) and b_triangular[i] != 0:
            Esol = False
            output_content.append(f"Fila ${i+1}$: sistema incompatible\n")
            display(Markdown(f"Fila ${i+1}$: sistema incompatible"))
            break
        elif all(A_triangular.row(i)[j] == 0 for j in range(m-1)) and b_triangular[i] == 0:
            inf_sol = True
            output_content.append(f"Fila ${i+1}$: \n\n ${True}$\n")
            display(Markdown(f"Fila ${i+1}$:\n\n ${latex(True)}$\n"))

    if not Esol:
        output_content.append("Sistema incompatible.\n")
        display(Markdown("Sistema incompatible."))
        return

    soluciones = {}
    variables = symbols(f'x1:{m}')
    columnas_cero = [j for j in range(m-1) if all(A_triangular[:, j][i] == 0 for i in range(n))]
    
    if columnas_cero_indices:
        # Identificar filas con no ceros y expresarlas como ecuaciones
        for i in range(n):
            if any(A_triangular.row(i)[j] != 0 for j in range(m-1)):
                ecuacion = Eq(A_triangular.row(i).dot(Matrix(variables[:m-1])), b_triangular[i])
                display(Markdown(f"Ecuación {i+1}:"))
                display(ecuacion)
                output_content.append(f"Ecuación {i+1}:\n{latex(ecuacion)}\n")

                # Resolver la ecuación
                solucion = solve(ecuacion, dict=True)
                if solucion:
                    # Despejar soluciones para variables dependientes
                    solucion = solucion[0]  # Tomar la primera solución si hay varias
                    for var, val in solucion.items():
                        if var not in soluciones:
                            soluciones[var] = val
                            display(Markdown(f"Despejamos para obtener:\n\n ${latex(var)} = {latex(val)}$"))
                            output_content.append(f"Solución para ${latex(var)}$: ${latex(val)}$\n")
                else:
                    display(Markdown(f"No se encontró solución para la ecuación {i+1}"))
                    output_content.append(f"No se encontró solución para la ecuación {i+1}\n")

    elif inf_sol:
        A_rref, b_rref = rrefsymb(matriz_ampliada)[:, :-1], rrefsymb(matriz_ampliada)[:, -1]
        display(Markdown("---"))

        for i in reversed(range(n)):
            ecuacion = Eq(A_rref.row(i).dot(Matrix(variables[:m-1])), b_rref[i])
            solucion = solve(ecuacion, variables[i])
            display(Markdown(f"Ecuación original ${i+1}$:"))
            display(ecuacion)
            output_content.append(f"Ecuación original ${i+1}$:\n{latex(ecuacion)}\n")
            if solucion:
                solucion_despejada = solucion[0] if isinstance(solucion, list) else solucion
                soluciones[variables[i]] = solucion_despejada

        # Agregar soluciones para variables libres (columnas de ceros)
        for indice in columnas_cero_indices:
            if variables[indice] not in soluciones:
                soluciones[variables[indice]] = "cualquier valor real"
                display(Markdown(f"${latex(variables[indice])}$ puede tomar cualquier valor real."))
                output_content.append(f"${latex(variables[indice])}$ puede tomar cualquier valor real.\n")
    else:
        for i in reversed(range(n)):
            ecuacion = Eq(A_triangular.row(i).dot(Matrix(variables[:m-1])), b_triangular[i])
            solucion = solve(ecuacion, variables[i])

            display(Markdown(f"Ecuación original {i+1}:"))
            display(ecuacion)
            output_content.append(f"Ecuación original {i+1}:\n{latex(ecuacion)}\n")      
            
            if solucion:
                solucion_despejada = solucion[0] if isinstance(solucion, list) else solucion
                display(Markdown(f"Despeje de $x_{i+1}$:"))
                display(Eq(variables[i], solucion_despejada))
                output_content.append(f"Despeje de $x_{i+1}$:\n{latex(Eq(variables[i], solucion_despejada))}\n")

                solucion_sustituida = solucion_despejada
                for var in soluciones:
                    solucion_sustituida = solucion_sustituida.subs(var, soluciones[var])
                soluciones[variables[i]] = solucion_sustituida
                display(Markdown(f"Solución sustituida para $x_{i+1}$:"))
                display(Eq(variables[i], solucion_sustituida))
                output_content.append(f"Solución sustituida para $x_{i+1}$:\n{latex(Eq(variables[i], solucion_sustituida))}\n")

    display(Markdown("Por lo tanto, las soluciones son:"))
    output_content.append("Por lo tanto, las soluciones son:\n")

    # Procesar primero las columnas con ceros para agregarlas a las soluciones
    for j in columnas_cero:
        variable = symbols(f"x_{j+1}")
        solucion_latex = latex(variable) + r" \in \mathbb{R}"  # Uso de r antes de la cadena para tratarlo como raw string
        soluciones[variable] = solucion_latex

    # Ordenar las soluciones por la variable
    soluciones_ordenadas = dict(sorted(soluciones.items(), key=lambda item: str(item[0])))

    # Imprimir y almacenar las soluciones
    for var, sol in soluciones_ordenadas.items():
        if isinstance(sol, str):
            # Directamente mostrar la solución si es una cadena (para columnas de ceros)
            display(Markdown(f"${latex(var)} = {sol}$"))
            output_content.append(f"${latex(var)} = {sol}$\n")
        else:
            # Convertir la solución a LaTeX en modo 'inline' para evitar el centrado y la numeración
            display(Markdown(f"${latex(var)} = {latex(sol)}$"))
            output_content.append(f"${latex(var)} = {latex(sol)}$\n")

    with open(archivo_salida, "w", encoding='utf-8') as file:
        file.writelines(output_content)

    return soluciones_ordenadas

# solve_GJ_alt: Descripción no disponible.
def solve_GJ_alt(A, b, archivo_salida="solucion_sistema.txt"):
    
    # Asumir que A ya es una matriz triangular inferior
    A_triangular = Matrix(A)
    b_triangular = Matrix(b)

    n, m = A_triangular.shape
    output_content = []  # Lista para almacenar el contenido que se imprime y guarda


    Esol = True
    inf_sol = False

    # Identificar columnas de ceros
    columnas_cero_indices = [j for j in range(m) if all(A_triangular[:, j][i] == 0 for i in range(n))]
    for i in range(n):
        if all(A_triangular.row(i)[j] == 0 for j in range(m-1)) and b_triangular[i] != 0:
            Esol = False
            output_content.append(f"Fila {i+1}: sistema incompatible\n")
            display(Markdown(f"Fila {i+1}: sistema incompatible"))
            break
        elif all(A_triangular.row(i)[j] == 0 for j in range(m-1)) and b_triangular[i] == 0:
            inf_sol = True
            output_content.append(f"Fila {i+1}: infinitas soluciones\n")
            display(Markdown(f"Fila {i+1}: infinitas soluciones"))

    if not Esol:
        output_content.append("Sistema incompatible.\n")
        display(Markdown("Sistema incompatible."))
        return

    soluciones = {}
    variables = symbols(f'y1:{m+1}')

    if columnas_cero_indices:
        for i in range(n):
            if any(A_triangular.row(i)[j] != 0 for j in range(m-1)):
                ecuacion = Eq(A_triangular.row(i).dot(Matrix(variables[:m])), b_triangular[i])
                display(Markdown(f"Ecuación {i+1}:"))
                display(ecuacion)
                output_content.append(f"Ecuación {i+1}:\n{latex(ecuacion)}\n")

                solucion = solve(ecuacion, dict=True)
                if solucion:
                    solucion = solucion[0]
                    for var, val in solucion.items():
                        if var not in soluciones:
                            soluciones[var] = val
                            display(Markdown(f"Despejamos para obtener:\n\n ${latex(var)} = {latex(val)}$"))
                            output_content.append(f"Solución para ${latex(var)}$: ${latex(val)}$\n")
                else:
                    display(Markdown(f"No se encontró solución para la ecuación {i+1}"))
                    output_content.append(f"No se encontró solución para la ecuación {i+1}\n")
        # Agregar soluciones para variables libres
        for indice in columnas_cero_indices:
            if variables[indice] not in soluciones:

                display(Markdown(f"${latex(variables[indice])}$ puede tomar cualquier valor real."))
                output_content.append(f"${latex(variables[indice])}$ puede tomar cualquier valor real.\n")
    else:
        for i, var in enumerate(variables[:n]):
            if var not in soluciones:
                # Sustituir soluciones anteriores en la ecuación actual
                ecuacion_actual = Eq(A_triangular.row(i).dot(Matrix(variables[:m])), b_triangular[i])
                for var_previa in soluciones:
                    ecuacion_actual = ecuacion_actual.subs(var_previa, soluciones[var_previa])

                display(Markdown(f"Solución sustituida para ${latex(var)}$:"))
                display(ecuacion_actual)
                output_content.append(f"Solución sustituida para ${latex(var)}$:\n{latex(ecuacion_actual)}\n")

                solucion = solve(ecuacion_actual, var)
                if solucion:
                    solucion_despejada = solucion[0] if isinstance(solucion, list) else solucion
                    soluciones[var] = solucion_despejada
                    display(Markdown(f"${latex(var)} = {latex(solucion_despejada)}$"))
                    output_content.append(f"${latex(var)} = {latex(solucion_despejada)}$\n")

    # Procesar primero las columnas con ceros para agregarlas a las soluciones
    for indice in columnas_cero_indices:
        variable = symbols(f"y_{indice+1}")
        solucion_latex = latex(variable) + r" \in \mathbb{R}"
        soluciones[variable] = solucion_latex

    # Ordenar las soluciones por la variable
    soluciones_ordenadas = dict(sorted(soluciones.items(), key=lambda item: str(item[0])))

    # Imprimir las soluciones
    display(Markdown("Por lo tanto, las soluciones son:"))
    output_content.append("Por lo tanto, las soluciones son:\n")

    # Imprimir y almacenar las soluciones
    for var, sol in soluciones_ordenadas.items():
        if isinstance(sol, str):
            display(Markdown(f"${latex(var)} = {sol}$"))
            output_content.append(f"${latex(var)} = {sol}$\n")
        else:
            display(Markdown(f"${latex(var)} = {latex(sol)}$"))
            output_content.append(f"${latex(var)} = {latex(sol)}$\n")

    with open(archivo_salida, "w", encoding='utf-8') as file:
        file.writelines(output_content)

    return soluciones_ordenadas

# comprobar_li: Descripción no disponible.
def comprobar_li(M, reduccion=0, archivo_salida="comprobar_li.txt"):

    A = Matrix(M)
    n, m = A.shape
    output_content = []  # Lista para almacenar el contenido impreso

    if reduccion == 0:
        # Usar rrefsymb para comprobar LI
        display(Markdown("Reducimos por filas :\n"))
        output_content.append("Reducimos por filas :\n\n")
        display(Markdown("---"))
        A_rref = rrefsymb(A)
        display(Markdown("---"))
        base = [A_rref.row(i) for i in range(n) if any(A_rref.row(i))]

    elif reduccion == 1 or reduccion == 2:
        # Usar reduccion_triangular para comprobar LI
        display(Markdown("Reducimos por filas ::\n"))
        output_content.append("Reducimos por filas ::\n\n")
        A_triangular = reduccion_triangular(A)
        display(Markdown("---"))
        if reduccion == 1:
            base = [A_triangular.row(i) for i in range(n) if any(A_triangular.row(i))]
        else:
            base = [A.row(i) for i in range(n) if any(A_triangular.row(i))]
        display(Markdown("---"))

    # Mostrar y guardar la base
    display(Markdown("Por lo que una base de las filas LI es:\n"))
    output_content.append("Por lo que una base de las filas LI es:\n\n")
    ind = 1
    for fila in base:
        display(Markdown(f"$v_{ind} = {latex(fila)}$"))
        output_content.append(f"$v_{ind} = {latex(fila)}$\n")
        ind += 1


    # Guardar en archivo
    with open(archivo_salida, "w", encoding='utf-8') as file:
        file.writelines(output_content)

    return base

# calcular_espacio_fila: Descripción no disponible.
def calcular_espacio_fila(M, archivo_salida="espacio_fila.txt"):

    A = Matrix(M)
    n, m = A.shape
    output_content = []  # Lista para almacenar el contenido impreso

    # Obtener la base del espacio fila
    base = comprobar_li(A, reduccion=0)
    display(Markdown("La base del espacio fila es:\n"))
    output_content.append("La base del espacio fila es:\n\n")
    for fila in base:
        display(Markdown(f"$v = {latex(fila)}$"))
        output_content.append(f"$v = {latex(fila)}$\n")

    # Expresar la base como una combinación lineal
    combinacion_total = ' + '.join([f"{latex(base[i].T)} x_{i + 1}" for i in range(len(base))])
    display(Markdown(f"Por lo tanto, el espacio fila se expresa como:\n\n $R(A) = {combinacion_total}$"))
    output_content.append(f"Por lo tanto, el espacio fila se expresa como:\n\n $R(A) = {combinacion_total}$\n")

    # Guardar en archivo
    with open(archivo_salida, "w", encoding='utf-8') as file:
        file.writelines(output_content)

    return base

# calcular_espacio_columna: Descripción no disponible.
def calcular_espacio_columna(M, archivo_salida="espacio_columna.txt"):

    A = Matrix(M).T  # Trabajar con la transpuesta de A
    n, m = A.shape
    output_content = []  # Lista para almacenar el contenido impreso

    # Obtener la base del espacio fila de A^T, que es la base del espacio columna de A
    base = comprobar_li(A, reduccion=0)
    display(Markdown("La base del espacio columna es:\n"))
    output_content.append("La base del espacio columna es:\n\n")
    for fila in base:
        display(Markdown(f"$v = {latex(fila)}$"))
        output_content.append(f"$v = {latex(fila)}$\n")

    # Expresar la base como una combinación lineal
    combinacion_total = ' + '.join([f"{latex(base[i].T)} x_{i + 1}" for i in range(len(base))])
    display(Markdown(f"Por lo tanto, el espacio columna se expresa como:\n\n $C(A) = {combinacion_total}$"))
    output_content.append(f"Por lo tanto, el espacio columna se expresa como:\n\n $C(A) = {combinacion_total}$\n")

    # Guardar en archivo
    with open(archivo_salida, "w", encoding='utf-8') as file:
        file.writelines(output_content)

    return base

# descomposicion_PLU: Descripción no disponible.
def descomposicion_PLU(A, archivo_salida="descomposicion_PLU.txt"):

    A = Matrix(A)
    n = A.rows
    U = A.copy()
    L = eye(n)
    P = eye(n)

    output_content = []  # Lista para almacenar todo el contenido que se imprime y guarda

    for i in range(n):
        # Encontrar el elemento máximo en la columna actual
        max_row = max(range(i, n), key=lambda x: abs(U[x, i]))
        if U[max_row, i] == 0:
            continue  # No hay elemento no nulo en esta columna

        # Intercambiar filas en U y P si es necesario
        if i != max_row:
            U.row_swap(i, max_row)
            P.row_swap(i, max_row)
            operacion = f"R_{{{i+1}}} \\leftrightarrow R_{{{max_row+1}}}"
            display(Markdown(f"Operación: ${operacion}$"))
            output_content.append(f"Operación: ${operacion}$\n")

        # Actualizar L con los multiplicadores de fila
        for j in range(i+1, n):
            L[j, i] = U[j, i] / U[i, i]
            U[j, :] = U[j, :] - L[j, i] * U[i, :]
            operacion = f"R_{{{j+1}}} - {L[j, i]}R_{{{i+1}}} \\leftarrow R_{{{j+1}}}"
            display(Markdown(f"Operación: ${operacion}$"))
            output_content.append(f"Operación: ${operacion}$\n")

        # Mostrar el paso actual
        display(Markdown(f"Paso {i+1}:"))
        display(Markdown("Matriz A:"))
        display(A)
        display(Markdown("Matriz L:"))
        display(L)
        display(Markdown("Matriz U:"))
        display(U)
        display(Markdown("Matriz P:"))
        display(P)

        # Agregar al contenido de salida
        output_content.append(f"Paso {i+1}:\n")
        output_content.append("Matriz A:\n" + latex(A) + "\n")
        output_content.append("Matriz L:\n" + latex(L) + "\n")
        output_content.append("Matriz U:\n" + latex(U) + "\n")
        output_content.append("Matriz P:\n" + latex(P) + "\n")

    # Guardar todo el contenido en un archivo en la carpeta raíz
    with open(archivo_salida, "w", encoding='utf-8') as file:
        file.writelines(output_content)

    return L, U, P

# solve_PLU: Descripción no disponible.
def solve_PLU(L, U, b, P=None, archivo_salida="solve_PLU.txt"):
    n = L.rows  # Asumiendo que L es cuadrada
    if P is None:
        P = eye(n)

    # Aplicar P a b
    Pb = Matrix(P) * Matrix(b)

    # Expresar la ecuación LUx = Pb
    display(Markdown(f"Resolviendo el sistema $LUx = Pb$ donde:"))
    display(Markdown(f"$L = {latex(L)}$, $U = {latex(U)}$, $P = {latex(P)}$, $b = {latex(Matrix(b))}$"))

    # Expresar Pb
    display(Markdown(f"$Pb = {latex(Pb)}$"))
    display(Markdown("---"))
    
    # Expresar el sistema Ly = Pb
    display(Markdown("Resolviendo sistema $Ly = Pb$:"))
    y = solve_GJ_alt(L, Pb, archivo_salida="solve_PLU_Ly.txt")
    display(Markdown("---"))
    
    # Transformar las soluciones en un arreglo
    y_sol = soluciones_a_arreglo(y)
    
    # Expresar el sistema Ux = y
    display(Markdown("Resolviendo sistema $Ux = y$:"))
    x = solve_GJ(U, y_sol, archivo_salida="solve_PLU_Ux.txt", triang=True)

    display(Markdown("---"))

    # Guardar las soluciones en un archivo
    with open(archivo_salida, "w", encoding='utf-8') as file:
        file.write("Soluciones obtenidas de L y U:\n")
        file.write(f"y = {latex(Matrix(y_sol))}\n")
        file.write(f"x = {latex(Matrix(soluciones_a_arreglo(x)))}\n")

    return x

# soluciones_a_arreglo: Descripción no disponible.
def soluciones_a_arreglo(soluciones):
    # Convertir las claves Symbol a cadenas y asegurarse de que estén en el orden correcto
    claves_ordenadas = sorted([str(clave) for clave in soluciones.keys()], key=lambda x: int(x[1:]))

    # Crear un arreglo con los valores en el orden correcto
    arreglo_soluciones = [soluciones[Symbol(clave)] for clave in claves_ordenadas]

    return arreglo_soluciones

# factorizar_polinomio: Descripción no disponible.
def factorizar_polinomio(A, lambda_symbol):
    # Polinomio característico y valores propios
    polinomio_caracteristico = A.charpoly(lambda_symbol).as_expr()
    valores_propios = solve(Eq(polinomio_caracteristico, 0), lambda_symbol)

    # Factorización utilizando los valores propios
    factores = [lambda_symbol - valor for valor in valores_propios]
    polinomio_factorizado = Mul(*factores)

    return polinomio_factorizado
# valores_propios: Descripción no disponible.
def valores_propios(M, archivo_salida="VALORES_PROPIOS.txt", lista = False):

    A = Matrix(M)
    lambda_symbol = symbols('lambda')
    n, m = A.shape

    output_content = []  # Lista para almacenar todo el contenido que se imprime y guarda

    if n != m:
        display(Markdown("La matriz debe ser cuadrada"))
        output_content.append("La matriz debe ser cuadrada\n")
        return

    # Polinomio característico y valores propios
    polinomio_caracteristico = A.charpoly(lambda_symbol)
    valores_propios = A.eigenvals()

    # Ecuación característica
    ecuacion_caracteristica = Eq(polinomio_caracteristico.as_expr(), 0)
    display(Markdown("**Ecuación característica:**"))
    display(ecuacion_caracteristica)
    output_content.append("**Ecuación característica:**\n")
    output_content.append(f"{latex(ecuacion_caracteristica)}\n\n")

    # Factorización con soluciones
    factorizacion = factorizar_polinomio(A, lambda_symbol)
    display(Markdown("**Factorización del polinomio característico:**"))
    display(factorizacion)
    output_content.append("**Factorización del polinomio característico:**\n")
    output_content.append(f"{latex(factorizacion)}\n\n")

    display(Markdown("---"))
    output_content.append("---\n")
    
    # Valores propios
    display(Markdown("**Valores propios:**"))
    output_content.append("**Valores propios:**\n")
    i = 1
    for valor, multiplicidad in valores_propios.items():
        display(Markdown(f"Valor propio: $\\lambda_{{{i}}} = {latex(valor)}$"))
        display(Markdown(f"Multiplicidad: {multiplicidad}"))
        output_content.append(f"Valor propio: $\\lambda_{{{i}}} = {latex(valor)}$\n")
        output_content.append(f"Multiplicidad: {latex(multiplicidad)}\n")
        output_content.append("---\n")
        i += 1

    # Guardar todo el contenido en un archivo en la carpeta raíz
    with open(archivo_salida, "w", encoding='utf-8') as file:
        file.writelines(output_content)
        
    
    if (lista): 
      valores_propios_lista = [valor for valor, _ in A.eigenvals().items()]
      return valores_propios_lista
    
    return valores_propios

# vectores_propios: Descripción no disponible.
def vectores_propios(M, archivo_salida="VECTORES_PROPIOS.txt", ant = True, uni = False):

    A = Matrix(M)
    lambda_symbol = symbols('lambda')
    n, m = A.shape

    output_content = []  # Lista para almacenar todo el contenido que se imprime y guarda
    
    if (ant): valores_propios(M, archivo_salida="VALORES_PROPIOS.txt",)

    if n != m:
        display(Markdown("La matriz debe ser cuadrada"))
        output_content.append("La matriz debe ser cuadrada\n")
        return

    # Calcula los valores propios y vectores propios
    valores_y_vectores_propios = A.eigenvects()
    vectores_propios_resultados = []

    for valor, multiplicidad, vectores in valores_y_vectores_propios:
        display(Markdown(f"Vectores propios asociados a $\\lambda = {latex(valor)}$:"))
        output_content.append(f"Vectores propios asociados a $\\lambda = {latex(valor)}$:\n")

        # Calcula A - lambda * I
        A_lambda_I = A - valor * Matrix.eye(n)
        display(Markdown(f"Matriz $A - \\lambda I = A - ({latex(valor)})I$:"))
        display(A_lambda_I)
        output_content.append("Matriz $A - \\lambda I$:\n")
        output_content.append(f"{latex(A_lambda_I)}\n")

        # Reducción de la matriz
        A_lambda_I_reducida, pivotes = A_lambda_I.rref()
        display(Markdown("Reducción de la matriz $A - \\lambda I$:"))
        display(A_lambda_I_reducida)
        output_content.append("Reducción de la matriz $A - \\lambda I$:\n")
        output_content.append(f"{latex(A_lambda_I_reducida)}\n")

        # Sistema de ecuaciones resultante
        display(Markdown("Sistema de ecuaciones resultante:"))
        output_content.append("Sistema de ecuaciones resultante:\n")
        variables = symbols('x1:%d' % (n+1))  # Crea variables x1, x2, ..., xn
        sistema_ecuaciones = [simplify(sum(row[i] * variables[i] for i in range(n))) for row in A_lambda_I_reducida.tolist()]
        ecuaciones_validas = [ec for ec in sistema_ecuaciones if ec != 0]
        for ecuacion in ecuaciones_validas:
            display(ecuacion)
            output_content.append(f"{latex(ecuacion)}\n")

        # Despeje de variable dependiente
        if ecuaciones_validas:
            soluciones = solve(ecuaciones_validas, variables, dict=True)
            display(Markdown("Despeje de variable dependiente:"))
            output_content.append("Despeje de variable dependiente:\n")
            for solucion in soluciones:
                for var, valor in solucion.items():
                    display(Eq(var, valor))
                    output_content.append(f"{latex(Eq(var, valor))}\n")

        # Vectores propios
        display(Markdown("Vectores propios:"))
        output_content.append("Vectores propios:\n")
        for vector in vectores:
            display(vector)
            output_content.append(f"{latex(vector)}\n")
        display(Markdown("---"))
        output_content.append("---\n")

        vectores_propios_resultados.extend(vectores)

    # Guardar todo el contenido en un archivo en la carpeta raíz
    with open(archivo_salida, "w", encoding='utf-8') as file:
        file.writelines(output_content)
    
    if (uni): 
        vectores_propios_normalizados = [vec.normalized() for _, _, vecs in A.eigenvects() for vec in vecs]
        return vectores_propios_normalizados

    return vectores_propios_resultados

# diago_vec_prop: Descripción no disponible.
def diago_vec_prop(M, archivo_salida="DIAGONALIZACION.txt", ant = False, inv = False):

    A = Matrix(M)
    lambda_symbol = symbols('lambda')
    n, m = A.shape

    output_content = []  # Lista para almacenar todo el contenido que se imprime y guarda
    if(ant): vectores_propios(M, archivo_salida="VECTORES_PROPIOS.txt", ant = True)
    if n != m:
        display(Markdown("La matriz debe ser cuadrada"))
        output_content.append("La matriz debe ser cuadrada\n")
        return

    # Calcula los valores propios y vectores propios
    valores_y_vectores_propios = A.eigenvects()

    # Verificar que los valores propios sean diferentes entre sí
    if len(valores_y_vectores_propios) != n:
        display(Markdown("No todos los valores propios son distintos. La matriz no es diagonalizable."))
        output_content.append("No todos los valores propios son distintos. La matriz no es diagonalizable.\n")
        return

    # Construir la matriz diagonal D con los valores propios
    D = Matrix.diag(*[val for val, mult, _ in valores_y_vectores_propios])
    display(Markdown("Matriz diagonal $D$:"))
    display(D)
    output_content.append("Matriz diagonal $D$:\n")
    output_content.append(f"{latex(D)}\n")

    # Construir la matriz X de vectores propios
    X = Matrix.hstack(*[vec[0] for _, _, vec in valores_y_vectores_propios])
    display(Markdown("Matriz de vectores propios $X$:"))
    display(X)
    output_content.append("Matriz de vectores propios $X$:\n")
    output_content.append(f"{latex(X)}\n")

    # Construir X^{-1}
    display(Markdown("---"))
    X_inv = X.inv()
    if(inv): inversa(X, archivo_salida="INVERSA_Diago.txt")
    display(Markdown("---"))
    display(Markdown("Matriz inversa de $X$, $X^{-1}$:"))
    display(X_inv)
    output_content.append("Matriz inversa de $X$, $X^{-1}$:\n")
    output_content.append(f"{latex(X_inv)}\n")

    # Mostrar la expresión XDX^{-1}
    display(Markdown("Expresión $XDX^{-1}$:"))
    display(Markdown(f"$XDX^{{-1}} = {latex(X)}{latex(D)}{latex(X_inv)}$"))
    output_content.append("Expresión $XDX^{-1}$:\n")
    output_content.append("Expresión $XDX^{-1}$: " + latex(X) + latex(D) + latex(X_inv) + "\n")

    # Guardar todo el contenido en un archivo en la carpeta raíz
    with open(archivo_salida, "w", encoding='utf-8') as file:
        file.writelines(output_content)

    return D, X
   
# printMsym: Descripción no disponible.
def printMsym(M):
    M = [Matrix(m) for m in M]
    combined_matrix = Matrix.hstack(*M)
    matrix_latex = latex(combined_matrix, mode='inline')
    display(Markdown(f"{matrix_latex}"))
    return matrix_latex

# gram_schmidt: Descripción no disponible.
def gram_schmidt(M, p_bas=False, p_nor=True, nombre_del_archivo = "gram_schmidt.txt"):    
    M = Matrix(M)
    n, m = M.shape
    U = []
    output_content = []  # Lista para almacenar todo el contenido que se imprime

    for i in range(m):
        u = M[:, i]
        for j in range(i):
            pp = U[j].dot(u)
            norma = U[j].dot(U[j])
            proj = (pp / norma) * U[j]
            u -= proj

            # Acumular el contenido a guardar
            output_content.append(f"Producto punto de $u_{j+1}$ y $v_{i+1}$: ${pp}$\n")
            output_content.append(f"Norma de $u_{j+1}$: ${norma}$\n")
            output_content.append(f"Proyección de v_{i+1} sobre $u_{j+1}$:\n${latex(proj)}$\n")
            output_content.append(f"$ Pr_{{u_{i}}} = {latex(proj)} $")
            output_content.append(f"Por lo que:\n\n $u_{i+1} = v_{i+1} - {latex(proj)} = {latex(u)}$")

            display(Markdown(f"Producto punto de $\\mathbf{{u}}_{{{j+1}}}$ y $\\mathbf{{v}}_{{{i+1}}}$: {pp}"))
            display(Markdown(f"Norma al cuadrado de $\\mathbf{{u}}_{{{j+1}}}$: {norma}"))
            display(Markdown(f"Proyección de $\\mathbf{{v}}_{{{i+1}}}$ sobre $\\mathbf{{u}}_{{{j+1}}}$:"))
            display(Markdown(f"$ Pr_{{u_{i}}} = {latex(proj)} $"))
            display(Markdown(f"Por lo que:\n\n $u_{i+1} = v_{i+1} - {latex(proj)} = {latex(u)}$"))

        U.append(u)
        if p_bas:
            output_content.append(f"Vector ortogonal $u_{i+1}$ antes de la normalización:\n")
            output_content.append(f"$u_{i+1} = {latex(u)}$")
            display(Markdown(f"Vector ortogonal $\\mathbf{{u}}_{{{i+1}}}$ antes de la normalización:"))
            display((Markdown(f"$u_{i+1} = {latex(u)}$")))

        output_content.append("---\n")
        display(Markdown("---"))

    display(Markdown("Base ortogonal:"))
    output_content.append("Base ortogonal:\n")
    orthogonal_base_latex = printMsym(U)
    output_content.append(orthogonal_base_latex + "\n")

    U_norm = []
    norms = []
    for u in U:
        norm = sqrt(u.dot(u))
        norms.append(norm)
        if norm != 0:
            U_norm.append(u / norm)
        else:
            U_norm.append(u)

    display(Markdown("Normas de los vectores ortogonales:"))
    output_content.append("Normas de los vectores ortogonales:\n")
    i = 1
    for norm in norms:
        output_content.append(f"Norma de $u_{i}$:\n")
        output_content.append(f"$||u_{i}|| = {latex(norm)}$")
        display(Markdown(f"Norma de $\\mathbf{{u}}_{{{i}}}$:"))
        display(Markdown(f"$||u_{i}|| = {latex(norm)}$"))
        i += 1

    display(Markdown("Base ortonormal:"))
    output_content.append("Base ortonormal:\n")
    orthonormal_base_latex = printMsym(U_norm)
    output_content.append(orthonormal_base_latex + "\n")

    # Guardar todo el contenido en un archivo
    # Guardar todo el contenido en un archivo
    with open(nombre_del_archivo, "w", encoding='utf-8') as file:
        file.writelines(output_content)


    return U_norm

# vectores_a_matriz_sympy: Descripción no disponible.
def vectores_a_matriz_sympy(vectores):
    # Convertir los vectores en una sola matriz de SymPy
    matriz = Matrix.hstack(*[Matrix(v) for v in vectores])
    return matriz

# factorizacion_QR: Descripción no disponible.
def factorizacion_QR(A, archivo_salida="factorizacion_QR.txt"):
    # Convertir A en una matriz de SymPy
    A = Matrix(A)
    display(Markdown("**Gram_Schmidt**:"))
    # Obtener la matriz Q usando Gram-Schmidt
    Q_list = gram_schmidt(A, p_bas=False, p_nor=True)
    Q = vectores_a_matriz_sympy(Q_list)

    display(Markdown("---"))
    # Calcular R como Q^T * A
    R = Q.T * A

    # Guardar la información en un archivo
    with open(archivo_salida, "w", encoding='utf-8') as file:
        file.write("Factorización QR de la matriz A:\n")
        file.write(f"Matriz $A:\\\\{latex(A)}\n\n")
        file.write("Matriz $Q$ (ortonormal):\n")
        file.write(f"{latex(Q)}\n\n")
        file.write("Matriz R (Q^T * A):\n")
        file.write(f"$R = Q^T \\cdot A = {latex(R)}$\n")

    # Mostrar la información
    display(Markdown("Factorización $QR$ de la matriz $A$:"))
    display(Markdown(f"$A = {latex(A)}$"))
    display(Markdown("Matriz $Q$ (ortonormal):"))
    display(Markdown(f"$Q = {latex(Q)}$"))
    display(Markdown("Matriz $Q^T$ (Q transpuesta):"))
    display(Markdown(f"$Q^T = {latex(Q.T)}$"))
    display(Markdown("Matriz $R$ ($Q^T \\cdot A$):"))
    display(Markdown(f"$R = Q^T \\cdot A = {latex(R)}$"))

    return Q, R

# minimos_cuadrados: Descripción no disponible.
def minimos_cuadrados(A, b, archivo_salida="minimos_cuadrados.txt"):
    # Convertir A y b en matrices de SymPy
    A = Matrix(A)
    b = Matrix(b)
    
    display(Markdown("Resolvemos por minimos cuadrados $Ax = b$:"))
    display(Markdown(f"${latex(A)} x = {latex(b)}$"))

    # Calcular y expresar A^T
    A_T = A.T
    display(Markdown("Matriz transpuesta $A^T$:"))
    display(A_T)
    
    display(Markdown("Buscamos resolver $A^T A = A^T b$:"))
    display(Markdown(f"${latex(A_T)} {latex(A)} x = {latex(A_T)} {latex(b)}$"))

    # Calcular y expresar A^T * A
    A_T_A = A_T * A
    display(Markdown("Matriz $A^T A$:"))
    display(A_T_A)
    

    # Calcular y expresar A^T * b
    A_T_b = A_T * b
    display(Markdown("Matriz $A^T b$:"))
    display(A_T_b)
    
    display(Markdown("Así $A^T A = A^T b$ es:"))
    display(Markdown(f"${latex(A_T_A)} x = {latex(A_T_b)}$"))
    
    # Resolver el sistema (A^T A)x = A^T b
    x = solve_GJ(A_T_A, A_T_b)

    # Guardar la información en un archivo
    with open(archivo_salida, "w", encoding='utf-8') as file:
        file.write("Resolvemos por minimos cuadrados $Ax = b$:\n\n")
        file.write(f"${latex(A)} x = {latex(b)}$\n\n")
        file.write(f"Matriz $A = {latex(A)}$\n\n")
        file.write(f"Matriz $b = {latex(b)}$\n\n")
        
        file.write("Buscamos resolver $A^T A = A^T b$:\n\n")     
        file.write(f"Matriz $A^T = {latex(A_T)}\n\n")
        file.write(f"Matriz $A^T A = {latex(A_T_A)}$\n\n")
        file.write(f"Matriz $A^T b = {latex(A_T_b)}$\n\n")
        file.write("Así $A^T A = A^T b$ es:\n\n")
        file.write(f"${latex(A_T_A)} x = {latex(A_T_b)}$")

        file.write(f"Solución de minimos cuadrados $x = {latex(x)}$\n")

    # Mostrar la solución
    display(Markdown("Solución $x$:"))
    display(x)

    return x

# minimos_cuadradosQR: Descripción no disponible.
def minimos_cuadradosQR(A, b, archivo_salida="minimos_cuadradosQR.txt"):
    # Convertir A y b en matrices de SymPy
    A = Matrix(A)
    b = Matrix(b)
    
    display(Markdown("Resolvemos por mínimos cuadrados $Ax = b$:"))
    display(Markdown(f"${latex(A)} x = {latex(b)}$"))

    # Factorización QR
    display(Markdown("---"))
    display(Markdown("Factorizacion $QR$:\n"))
    Q, R = factorizacion_QR(A, "gram_schmidtQR.txt")
    display(Markdown("---"))
    
    display(Markdown("Transponemos $Q^T = Q^{-1}$:\n"))
    display(Markdown(f"$Q^T = {latex(Q.T)}$\n"))
    
    # Calcular Q^T * b
    Q_T_b = Q.T * b
    
    display(Markdown("De forma que:\n"))
    display(Markdown(f"$Q^T b = {latex(Q_T_b)}$\n"))
    
    display(Markdown("Así resolvemos $Rx = Q^T b$:\n"))
    display(Markdown(f"${latex(R)} x = {latex(Q_T_b)}$\n"))

    # Resolver el sistema Rx = Q^Tb
    x = solve_GJ(R, Q_T_b, "sistema_sol_QR.txt", triang=True)

    # Guardar la información en un archivo
    with open(archivo_salida, "w", encoding='utf-8') as file:
        file.write("Resolvemos por mínimos cuadrados $Ax = b$:\n\n")
        file.write(f"${latex(A)} x = {latex(b)}$\n\n")
        file.write("---\n")
        file.write("Factorización QR de A:\n\n")
        file.write(f"Matriz $Q = {latex(Q)}$\n\n")
        file.write(f"Matriz $R = {latex(R)}$\n\n")
        file.write("---\n")
        file.write("Transponemos $Q^T = Q^{-1}$:\n\n")
        file.write(f"$Q^T = {latex(Q.T)}$\n\n")
        file.write("De forma que:\n\n")
        file.write(f"$Q^T b = {latex(Q_T_b)}$\n\n")
        file.write("Así resolvemos $Rx = Q^T b$:\n\n")
        file.write(f"${latex(R)} x = {latex(Q_T_b)}$\n\n")
        file.write("---\n")
        file.write("Solución de mínimos cuadrados $x = {latex(x)}$\n")

    # Mostrar la solución
    display(Markdown("Solución $x$:"))
    display(x)

    return x

# diago_vec_prop: Descripción no disponible.
def diago_vec_prop(M, archivo_salida="DIAGONALIZACION.txt", ant=False, inv=False):
    A = Matrix(M)
    lambda_symbol = symbols('lambda')
    n, m = A.shape

    output_content = []

    if ant: vectores_propios(M, archivo_salida="VECTORES_PROPIOS.txt", ant=True)

    if n != m:
        display(Markdown("La matriz debe ser cuadrada"))
        output_content.append("La matriz debe ser cuadrada\n")
        return

    valores_y_vectores_propios = A.eigenvects()

    # Calcular el número total de vectores propios
    total_vectores = sum(multiplicidad for _, multiplicidad, _ in valores_y_vectores_propios)

    if total_vectores < n:
        display(Markdown("No hay suficientes vectores propios para la diagonalización."))
        output_content.append("No hay suficientes vectores propios para la diagonalización.\n")
        return

    valores_y_vectores_propios = A.eigenvects()

    if len(valores_y_vectores_propios) != n:
        display(Markdown("No todos los valores propios son distintos. La matriz no es diagonalizable."))
        output_content.append("No todos los valores propios son distintos. La matriz no es diagonalizable.\n")
        return

    D = Matrix.diag(*[val for val, mult, _ in valores_y_vectores_propios])
    display(Markdown("Matriz diagonal $D$:"))
    display(D)
    output_content.append("Matriz diagonal $D$:\n")
    output_content.append(f"{latex(D)}\n")

    vectores_norm = [vec[0].normalized() for _, _, vec in valores_y_vectores_propios]
    P = Matrix.hstack(*vectores_norm)
    display(Markdown("Matriz de vectores propios $P$ normalizados:"))
    display(P)
    output_content.append("Matriz de vectores propios $P$ normalizados:\n")
    output_content.append(f"{latex(P)}\n")

    display(Markdown("---"))
    P_inv = P.inv()
    if inv: inversa(P, archivo_salida="INVERSA_Diago.txt")
    display(Markdown("---"))
    display(Markdown("Matriz inversa de $P$, $P^{-1}$:"))
    display(P_inv)
    output_content.append("Matriz inversa de $P$, $P^{-1}$:\n")
    output_content.append(f"{latex(P_inv)}\n")

    display(Markdown("Expresión $PDP^{-1}$:"))
    display(Markdown(f"$PDP^{{-1}} = {latex(P)}{latex(D)}{latex(P_inv)}$"))
    output_content.append("Expresión $PDP^{-1}$:\n")
    output_content.append("Expresión $PDP^{-1}$: " + latex(P) + latex(D) + latex(P_inv) + "\n")

    with open(archivo_salida, "w", encoding='utf-8') as file:
        file.writelines(output_content)

    # Regresar P, D y P^T
    return P, D

# forma_cuadratica_a_polinomio: Descripción no disponible.
def forma_cuadratica_a_polinomio(A, archivo_salida="forma_cuadratica_polinomio.txt"):
    n = A.shape[0]
    x = Matrix(symbols(f'x1:{n+1}'))  # Genera un vector de variables simbólicas x1, x2, ..., xn

    output_content = []  # Lista para almacenar todo el contenido que se imprime y guarda

    # Verificar si la matriz es simétrica
    if A != A.T:
        display(Markdown("La matriz debe ser simétrica"))
        output_content.append("La matriz debe ser simétrica\n")
        return

    # Expresar la forma cuadrática como x^T A x
    forma_cuadratica_expr = f"x^T {latex(A)} x"
    display(Markdown("La forma cuadrática es:"))
    display(Markdown(f"$Q(\\mathbf{{x}}) = {forma_cuadratica_expr}$"))
    output_content.append("La forma cuadrática es:\n")
    output_content.append(f"$Q(\\mathbf{{x}}) = {forma_cuadratica_expr}$\n")

    # Calcular la forma cuadrática
    forma_cuadratica = x.T * A * x

    # Expandir a forma polinomial
    polinomio = forma_cuadratica[0].expand()  # Acceder al primer elemento de la matriz resultante
    display(Markdown("Su forma polinomial es:"))
    display(Markdown(f"$Q(\\mathbf{{x}}) = {latex(polinomio)}$"))
    output_content.append("Su forma polinomial es:\n")
    output_content.append(f"$Q(\\mathbf{{x}}) = {latex(polinomio)}$\n")

    # Guardar todo el contenido en un archivo
    with open(archivo_salida, "w", encoding='utf-8') as file:
        file.writelines(output_content)

    return polinomio

# polinomio_a_forma_cuadratica: Descripción no disponible.
def polinomio_a_forma_cuadratica(polinomio_str, archivo_salida="polinomio_a_forma_cuadratica.txt"):
    # Extraer las variables del polinomio
    variables = sorted(list({str(s) for s in parse_expr(polinomio_str).free_symbols}), key=lambda v: v[1:])
    x = symbols(' '.join(variables))

    # Parsear el polinomio
    polinomio = parse_expr(polinomio_str, local_dict={var: x[i] for i, var in enumerate(variables)})

    # Inicializar una matriz de ceros
    n = len(x)
    A = Matrix.zeros(n)

    # Descomponer el polinomio y llenar la matriz A
    coef_dict = polinomio.as_coefficients_dict()
    for monomio, coef in coef_dict.items():
        indices = [i for i, xi in enumerate(x) if xi in monomio.free_symbols]

        # Llenar la matriz A
        if len(indices) == 2:
            i, j = indices
            A[i, j] = coef / 2  # Dividir por 2 los coeficientes no diagonales
            A[j, i] = coef / 2
        elif len(indices) == 1:
            i = indices[0]
            A[i, i] = coef

    # Expresar el polinomio
    display(Markdown("Su forma polinomial es:"))
    display(Markdown(f"$Q(\\mathbf{{x}}) = {latex(polinomio)}$"))

    # Expresar la forma cuadrática
    forma_cuadratica_expr = f"x^T {latex(A)} x"
    display(Markdown("La forma cuadrática es:"))
    display(Markdown(f"$Q(\\mathbf{{x}}) = {forma_cuadratica_expr}$"))

    # Guardar en un archivo
    with open(archivo_salida, "w", encoding='utf-8') as file:
        file.write("Su forma polinomial es:\n")
        file.write(f"$Q(\\mathbf{{x}}) = {latex(polinomio)}$\n\n")
        file.write("La forma cuadrática es:\n")
        file.write(f"$Q(\\mathbf{{x}}) = {forma_cuadratica_expr}$\n")

    return A

# verificar_definicion_matriz: Descripción no disponible.
def verificar_definicion_matriz(M, archivo_salida="verificar_definicion.txt"):
    # Obtener los valores propios de la matriz
    valores = valores_propios(M)

    # Inicializar variables para la conclusión
    conclusion = ""
    valores_positivos = [valor for valor in valores if valor > 0]
    valores_negativos = [valor for valor in valores if valor < 0]
    valores_cero = [valor for valor in valores if valor == 0]

    # Verificar la definición de la matriz
    todos_positivos = all(valor > 0 for valor in valores)
    todos_no_negativos = all(valor >= 0 for valor in valores)
    todos_negativos = all(valor < 0 for valor in valores)
    todos_no_positivos = all(valor <= 0 for valor in valores)

    # Preparar las partes de la cadena fuera de la f-string
    partes_lambda_positivos = [f"$\\lambda_{{{i+1}}} = {valor}$" for i, valor in enumerate(valores_positivos)]
    partes_lambda_negativos = [f"$\\lambda_{{{i+1}}} = {valor}$" for i, valor in enumerate(valores_negativos)]
    partes_lambda_cero = [f"$\\lambda_{{{i+1}}} = {valor}$" for i, valor in enumerate(valores_cero)]

    # Determinar la definición de la matriz y construir la conclusión
    if todos_positivos:
        definicion = "Positiva definida"
        conclusion = f"Todos los valores propios ({', '.join(partes_lambda_positivos)}) son mayores que $0$."
    elif todos_no_negativos:
        definicion = "No negativa definida"
        conclusion = f"Todos los valores propios ({', '.join(partes_lambda_cero)}) son no negativos."
    elif todos_negativos:
        definicion = "Negativa definida"
        conclusion = f"Todos los valores propios ({', '.join(partes_lambda_negativos)}) son menores que $0$."
    elif todos_no_positivos:
        definicion = "No positiva definida"
        conclusion = f"Todos los valores propios ({', '.join(partes_lambda_cero)}) son no positivos."
    else:
        definicion = "Indefinida"
        conclusion = f"La matriz tiene valores propios tanto positivos ({', '.join(partes_lambda_positivos)}) como negativos ({', '.join(partes_lambda_negativos)})."

    # Mostrar el resultado y la conclusión
    display(Markdown(f"La matriz es **{definicion}**.\n\n**Conclusión**: {conclusion}"))

    # Guardar en un archivo
    with open(archivo_salida, "w", encoding='utf-8') as file:
        file.write(f"La matriz es **{definicion}**.\n\n**Conclusión**: {conclusion}\n")

    return definicion

def construir_V(M, archivo_salida="Construir_V.txt", rref = False):

    A = Matrix(M)
    n, m = A.shape
    output_content = []  # Lista para almacenar el contenido impreso y guardado

    # Procesar A^TA
    display(Markdown("---"))
    display(Markdown("**Valores y vectores propios de $A^TA$**:\n"))
    ATA = A.transpose() * A
    display(Markdown(f"$A^TA = {latex(ATA)}$"))
    valores_propios(ATA, archivo_salida="VALORES_PROPIOS_ATA.txt")
    valores_y_vectores_propios_ATA = ATA.eigenvects()
    valores_y_vectores_propios_ATA.sort(key=lambda x: x[0], reverse=True)  # Ordena de mayor a menor
    vectores_propios_con_multiplicidad = []

    # Imprimir y almacenar los valores y vectores propios de A^TA
    ind = 1
    ind2 = 1
    for valor, multiplicidad, vectores in valores_y_vectores_propios_ATA:
        display(Markdown(f"**Valor propio:** $\\lambda = {latex(valor)}$"))
        output_content.append(f"Valor propio: \\lambda = ${latex(valor)}\n")
        
        A_lambda_I = ATA - valor * Matrix.eye(m)  # Matriz A^TA - lambda * I
        display(Markdown(f"Matriz $A^TA - \\lambda I = A^TA - ({latex(valor)})I$:"))
        display(Markdown(f"$A^TA - ({latex(valor)})I = {latex(A_lambda_I)}$"))
        output_content.append("Matriz $A^TA - \\lambda I$:\n")
        output_content.append(f"$A^TA - ({latex(valor)})I = {latex(A_lambda_I)}$")

        # Reducción de la matriz A^TA - lambda * I
        A_lambda_I_reducida, pivotes = A_lambda_I.rref()
        display(Markdown("Reducción de la matriz $A^TA - \\lambda I$:"))
        if rref: 
          display(Markdown("---"))
          rrefsymb(A_lambda_I, archivo_salida="RREF_SVD_V.txt")
          display(Markdown("---"))
        display(A_lambda_I_reducida)
        output_content.append("Reducción de la matriz $A^TA - \\lambda I$:\n")
        output_content.append(f"{latex(A_lambda_I_reducida)}\n")

        # Sistema de ecuaciones resultante
        display(Markdown("Sistema de ecuaciones resultante:"))
        output_content.append("Sistema de ecuaciones resultante:\n")
        variables = symbols('x1:%d' % (m+1))  # Crea variables x1, x2, ..., xm
        sistema_ecuaciones = [simplify(sum(row[i] * variables[i] for i in range(m))) for row in A_lambda_I_reducida.tolist()]
        ecuaciones_validas = [ec for ec in sistema_ecuaciones if ec != 0]
        for ecuacion in ecuaciones_validas:
            display(Markdown(f"${latex(ecuacion)} = 0$"))
            output_content.append(f"{latex(ecuacion)} = 0\n")

        # Despeje de variable dependiente
        if ecuaciones_validas:
            soluciones = solve(ecuaciones_validas, variables, dict=True)
            display(Markdown("Despeje de variable dependiente:"))
            output_content.append("Despeje de variable dependiente:\n")
            for solucion in soluciones:
                for var, valor in solucion.items():
                    display(Eq(var, valor))
                    output_content.append(f"{latex(Eq(var, valor))}\n")

        # Normalizar y almacenar TODOS los vectores propios asociados al valor propio
        for vector in vectores:
            vector_normalizado = vector.normalized()
            display(Markdown(f"**Vector propio de $A^TA$:** ${latex(vector)}$"))
            output_content.append(f"Vector propio de $A^TA$: ${vector_a_latex(vector)}$\n")
            display(Markdown(f"**Vector propio normalizado de $A^TA$:**"))
            display(Markdown(f"$v_{ind} = {latex(vector_normalizado)}$"))
            output_content.append(f"Vector propio normalizado de $A^TA$:\n")
            output_content.append(f"$v_{ind} = {vector_a_latex(vector_normalizado)}$\n")
            ind += 1

    # Formar matriz V
    vectores_propios_ATA_normalizados = [vector.normalized() for _, _, vectores in valores_y_vectores_propios_ATA for vector in vectores]
    V = Matrix.hstack(*vectores_propios_ATA_normalizados)
    display(Markdown(f"**Por lo tanto:**\n\n$V = {latex(V)}$"))
    output_content.append(f"Por lo tanto:\n\n$V = {latex(V)}$\n")
    
    # Guardar en archivo
    with open(archivo_salida, "w", encoding='utf-8') as file:
        file.writelines(output_content)
        
    return V, [sqrt(valor) for valor, _, _ in valores_y_vectores_propios_ATA if valor > 0]

def vectores_UNULL_numerico(AT_reducida, pivotes):
    # Dimensiones de la matriz transpuesta reducida
    n, m = AT_reducida.shape
    todas_las_variables = symbols('x1:%d' % (m+1))

    # Identificar variables libres (no pivotes)
    variables_libres = [i for i in range(m) if i not in pivotes]

    # Formar vectores del espacio nulo
    espacio_nulo_vectores = []
    for var_libre in variables_libres:
        vector_en_espacio_nulo = [0] * m  # Inicializar el vector con ceros
        vector_en_espacio_nulo[var_libre] = 1  # Asignar 1 a la variable libre actual

        # Resolver para las variables pivote
        for i in pivotes:
            suma = sum(AT_reducida[i, j] * todas_las_variables[j] for j in range(m))
            # Realizar sustituciones numéricas
            valor = -suma.subs({todas_las_variables[j]: (1 if j == var_libre else 0) for j in range(m)})
            vector_en_espacio_nulo[i] = valor

        vector_en_espacio_nulo = Matrix([vector_en_espacio_nulo]).T
        espacio_nulo_vectores.append(vector_en_espacio_nulo)

    return espacio_nulo_vectores

def gram_smith_SVD(Nulidad, r, output_content):
    U_ortogonalizados = []
    output_content.append(f"Aplicamos Gram Smith a los vectores del espacio nulo:\n")
    display(Markdown("Aplicamos Gram Smith a los vectores del espacio nulo:\n"))
    
    for i, a in enumerate(Nulidad):
        u = a
        for j in range(len(U_ortogonalizados)):
            u_j = U_ortogonalizados[j]
            pp = u.dot(u_j)
            norma = u_j.dot(u_j)
            proj = pp / norma * u_j
            u = u - proj

            # Acumular el contenido a guardar y mostrar
            output_content.append(f"Producto punto de $u_{j+1}$ y $a_{j + 2}$\n\n $u_{j+1} \\cdot a_{r+1+i} = {latex(pp)}$\n")
            output_content.append(f"Proyección de $a_{{{1 + i}}}$ sobre $u_{{{r + j + 1}}}$:")
            output_content.append(f"$Pr_{{u_{r + j + 1}}} a_{{{1 + i}}} = {vector_a_latex(proj)}$")
            display(Markdown(f"Producto punto de $u_{r + j + 1}$ y $a_{i + 2}$\n\n $u_{r + j + 1} \\cdot a_{1+i} = {latex(pp)}$\n"))
            display(Markdown(f"Proyección de $a_{{{1 + i}}}$ sobre $u_{{{r + j + 1}}}$:"))
            display(Markdown(f"$Pr_{{u_{r + j + 1}}} a_{{{1 + i}}} = {latex(proj)}$"))
        
        display(Markdown(f"Por lo que el vector ortogonal es: \n\n ${latex(u)}$"))
        output_content.append(f"Por lo que el vector ortogonal es: \n\n ${vector_a_latex(u)}$\n\n")
        
        norma_u = u.norm()
        if norma_u != 0:
            display(Markdown(f"Tenemos que: \n\n $||a_{i + 1}|| = {latex(norma_u)}$"))
            output_content.append(f"Tenemos que: \n\n $||a_{i + 1}|| = {latex(norma_u)}$")
            display(Markdown(f"Por lo que nomalizando obtenemos: \n\n $ u_{r+1+i} = {latex(u / norma_u)}$"))
            output_content.append(f"Por lo que: \n\n $u_{r+1+i} = {vector_a_latex(u / norma_u)}$\n\n")
            U_ortogonalizados.append(u / norma_u)

    display(Markdown("Por lo que los vectores faltantes son:\n"))
    output_content.append("Por lo que los vectores faltantes son:\n")
    
    for j, u in enumerate(U_ortogonalizados, start=r + 1):
        output_content.append(f"$u_{j} = {latex(u)}$\n")
        display(Markdown(f"$u_{j} = {latex(u)}$"))

    return Matrix.hstack(*[Matrix(v) for v in U_ortogonalizados])

def construir_U(A, V, Valores_Singulares, archivo_salida="Construir_U.txt", rref=False):
    n, m = A.shape
    output_content = []  # Lista para almacenar el contenido impreso y guardado

    display(Markdown("**Calculamos las columnas de $U$ usando $\\sigma_i u_i = Av_i$:**"))
    output_content.append("Calculamos las columnas de $U$ usando $\\sigma_i u_i = Av_i$:\n\n")

    
    # Calcular los u_i
    U_columnas = []
    for i in range(len(Valores_Singulares)):
        v_i = V.col(i)
        u_i = A * v_i / Valores_Singulares[i]
        U_columnas.append(u_i)
        display(Markdown(f"$u_{i + 1} = {latex(1 / Valores_Singulares[i])}Av_{i + 1}$\n\n"))
        display(Markdown(f"$u_{i + 1} = {latex(1 / Valores_Singulares[i])} {latex(A)} {latex(v_i)} = {latex(1 / Valores_Singulares[i])} {latex(A*v_i)} = {latex(u_i)} $\n\n"))
        output_content.append(f"$u_{i + 1} = {latex(1 / Valores_Singulares[i])} Av_{i + 1}$\n\n")
        output_content.append(f"$u_{i + 1} = {latex(1 / Valores_Singulares[i])} {latex(A)} {vector_a_latex(v_i)} = {latex(1 / Valores_Singulares[i])} {vector_a_latex(A*v_i)} = {vector_a_latex(u_i)} $\n\n")

        output_content.append(f" Por lo que $u_{i + 1} = {vector_a_latex(u_i)}$\n")
        display(Markdown(f" Por lo que $u_{i + 1} = {latex(u_i)}$\n"))


    # Completar U si r < n
    Nulidad = []
    if len(Valores_Singulares) < n:
        # Reducir la matriz A^T
        display(Markdown("Reducción de la matriz $A^T$:"))
        AT_reducida, pivotes = A.transpose().rref()
        output_content.append(f"Reducción de la matriz $A^T$\n\n $A^T = {latex(A.T)} \\rightarrow {latex(AT_reducida)}$\n")
        display(Markdown(f"$A^T = {latex(A.T)} \\rightarrow {latex(AT_reducida)}$"))

        # Sistema de ecuaciones resultante para el espacio nulo de A^T
        display(Markdown("Sistema de ecuaciones resultante para el espacio nulo de $A^T$:"))
        todas_las_variables = symbols('x1:%d' % (n+1))
        sistema_ecuaciones = [sum(row[i] * todas_las_variables[i] for i in range(n)) for row in AT_reducida.tolist()]
        ecuaciones_validas = [ec for ec in sistema_ecuaciones if ec != 0]
        for ecuacion in ecuaciones_validas:
            display(Markdown(f"${latex(ecuacion)} = 0$"))
            output_content.append(f"{latex(ecuacion)} = 0\n")

        # Despejar las variables dependientes
                # Despeje de variable dependiente
        if ecuaciones_validas:
            soluciones = solve(ecuaciones_validas, todas_las_variables, dict=True)
            display(Markdown("Despeje de variable dependiente:"))
            output_content.append("Despeje de variable dependiente:\n")
            for solucion in soluciones:
                for var, valor in solucion.items():
                    display(Eq(var, valor))
                    output_content.append(f"{latex(Eq(var, valor))}\n")
        
        num_vectores = n - len(Valores_Singulares)
        espacio_nulo_vectores = vectores_UNULL_numerico(AT_reducida, pivotes)
        k = 0
        display(Markdown("De esta forma tomamos los vectores del espacio nulo de $A^T$:"))
        output_content.append("De esta forma tomamos los vectores del espacio nulo de $A^T$:\n")

        for vector_en_espacio_nulo in espacio_nulo_vectores:
            Nulidad.append(vector_en_espacio_nulo)
            display(Markdown(f"$a_{k + 1} = {latex(vector_en_espacio_nulo)} $"))
            output_content.append(f"$a_{k + 1} = {vector_a_latex(vector_en_espacio_nulo)} $\n")

            if k + num_vectores == n:
                break
            k += 1

    # Gram Smith a U_columnas
    if Nulidad: 
        U_columnas.append(gram_smith_SVD(Nulidad, len(Valores_Singulares), output_content))
    
    # Formar matriz U con los vectores calculados
    U = Matrix.hstack(*U_columnas)
    display(Markdown("**De esta forma $U$ es**:\n\n"))
    display(Markdown(f"$U = {latex(U)}$"))
    output_content.append("De esta forma $U$ es:\n\n")
    output_content.append(f"$U = {latex(U)}$")
    # Guardar en archivo
    with open(archivo_salida, "w", encoding='utf-8') as file:
        file.writelines(output_content)

    return U

def svd_descomposicion(A, archivo_salida="svd_descomposicion.txt", rref=False):
    A = Matrix(A)
    n, m = A.shape
    output_content = []  # Lista para almacenar el contenido impreso

    display(Markdown("Haremos la descomposición SVD de la matriz:\n"))
    output_content.append("Haremos la descomposición SVD de la matriz:\n\n")
    display(Markdown(f"$A = {latex(A)}$\n"))
    output_content.append(f"$A = {latex(A)}$\n\n")

    display(Markdown("Buscaremos los vectores propios de:\n"))
    output_content.append("Buscaremos los vectores propios de:\n\n")
    AtA = A.T * A
    display(Markdown(f"$A^T A = {latex(AtA)}$\n"))
    output_content.append(f"$A^T A = {latex(AtA)}$\n")

    display(Markdown("---"))
    output_content.append("---\n")
    V, sing = construir_V(A, archivo_salida="Construir_V.txt", rref = False)
    display(Markdown("---"))
    output_content.append("---\n")
    
    display(Markdown("Los valores singulares son:\n"))
    output_content.append("Los valores singulares son:\n\n")
    i = 1
    for sigma in sing:
        display(Markdown(f"$\\sigma_{{{i}}} = {{{latex(sigma)}}}$\n"))
        output_content.append(f"$\\sigma_{{{i}}} = {{{latex(sigma)}}}$\n")
        i += 1

    # Crear D tal que sea diagonal y con los valores singulares distintos de cero
    D = diag(*sing)
    display(Markdown(f"$D = {latex(D)}$"))
    output_content.append(f"$D = {latex(D)}$\n")

    # Crear Sigma tal que sea nxm, poner D en la diagonal y puros ceros
    Sigma = zeros(n, m)
    for i in range(len(sing)):
        Sigma[i, i] = sing[i]
    display(Markdown(f"$\\Sigma = {latex(Sigma)}$"))
    output_content.append(f"$\\Sigma = {latex(Sigma)}$\n")

    display(Markdown("---"))
    output_content.append("---\n")
    U = construir_U(A, V, sing, archivo_salida="Construir_U.txt", rref=False)
    display(Markdown("---"))
    output_content.append("---\n")
    
    # Expresar la descomposicion SVD
    display(Markdown("**Por lo que la descomposición $SVD$ es**:\n"))
    output_content.append("Por lo que la descomposición SVD es:\n\n")
    display(Markdown(f"$A = U\\Sigma V^T = {latex(U)}{latex(Sigma)}{latex(V.T)}$"))
    output_content.append(f"$A = U\\Sigma V^T = {latex(U)}{latex(Sigma)}{latex(V.T)}$\n")

    # Guardar en archivo
    with open(archivo_salida, "w", encoding='utf-8') as file:
        file.writelines(output_content)

    return U, Sigma, V.T

# svd_reducida: Descripción no disponible.
def svd_reducida(U, S, VT, archivo_salida="svd_reducida.txt"):
    n, m = S.shape
    output_content = []  # Lista para almacenar el contenido impreso
    
    display(Markdown("---"))
    output_content.append("---\n")
    display(Markdown("**Forma Reducida**"))
    output_content.append("Forma Reducida\n")
    display(Markdown("---"))
    output_content.append("---\n")
     
    # Determinar el rango de la matriz S (número de valores singulares no nulos)
    rango = sum([1 for valor in S.diagonal() if valor > 0])
    display(Markdown(f"El rango de la matriz $\\Sigma$ es {rango}.\n"))
    output_content.append(f"El rango de la matriz $S$ es {rango}.\n\n")

    # Reducir U, S y VT al rango r
    U_reducida = U[:, :rango]
    S_reducida = S[:rango, :rango]
    VT_reducida = VT[:rango, :]

    display(Markdown(f"La matriz $U$ reducida es:\n\n$U_r = {latex(U_reducida)}$"))
    output_content.append(f"La matriz $U$ reducida es:\n\n$U_r = {latex(U_reducida)}$\n\n")
    display(Markdown(f"La matriz $S$ reducida es:\n\n$\\Sigma_r = {latex(S_reducida)}$"))
    output_content.append(f"La matriz $S$ reducida es:\n\n$\\Sigma_r = {latex(S_reducida)}$\n\n")
    display(Markdown(f"La matriz $V^T$ reducida es:\n\n$V_r^T = {latex(VT_reducida)}$"))
    output_content.append(f"La matriz $V^T$ reducida es:\n\n$V_r^T = {latex(VT_reducida)}$\n\n")
    
    

    # Guardar en archivo
    with open(archivo_salida, "w", encoding='utf-8') as file:
        file.writelines(output_content)

    return U_reducida, S_reducida, VT_reducida

# inversa_moore_penrose: Descripción no disponible.
def inversa_moore_penrose(U_r, S_r, VT_r, archivo_salida="inversa_moore_penrose.txt"):
    
    U_r = Matrix(U_r)
    S_r = Matrix(S_r)
    VT_r = Matrix(VT_r)

    output_content = []  # Lista para almacenar el contenido impreso

    # V_r es la transpuesta de VT_r
    V_r = VT_r.T
    display(Markdown(f"Como $V^T_r = {latex(VT_r)}$ entonces"))
    display(Markdown(f"$V_r = {latex(V_r)}$"))
    output_content.append(f"Como $V^T_r = {latex(VT_r)}$ entonces\n")
    output_content.append(f"$V_r = {latex(V_r)}$\n\n")

    # Calcular la inversa de S_r
    S_r_inv = diag(*[1/s if s != 0 else 0 for s in S_r.diagonal()])
    display(Markdown(f"Como $S_r = {latex(S_r)}$ entonces"))
    display(Markdown(f"$S_r^{{-1}} = {latex(S_r_inv)}$"))
    output_content.append(f"Como $S_r = {latex(S_r)}$ entonces\n")
    output_content.append(f"$S_r^{{-1}} = {latex(S_r_inv)}$\n\n")

    # U_r^T es la transpuesta de U_r
    U_r_T = U_r.T
    display(Markdown(f"Como $U_r = {latex(U_r)}$ entonces"))
    display(Markdown(f"$U_r^T = {latex(U_r_T)}$"))
    output_content.append(f"Como $U_r = {latex(U_r)}$ entonces\n")
    output_content.append(f"$U_r^T = {latex(U_r_T)}$\n\n")

    # Calcular la inversa de Moore-Penrose
    A_dag = V_r * S_r_inv * U_r_T
    display(Markdown("Por lo que la inversa de Moore-Penrose $A^{\\dag}$ es:"))
    display(Markdown(f"$A^{{\\dag}} = {latex(V_r)} {latex(S_r_inv)} {latex(U_r_T)}$"))
    display(Markdown(f"$A^{{\\dag}} = {latex(V_r * S_r_inv)} {latex(U_r_T)}$"))
    display(Markdown(f"$A^{{\\dag}} = {latex(A_dag)}$"))
    output_content.append("Por lo que la inversa de Moore-Penrose $A^{\\dag}$ es:\n")
    output_content.append(f"$A^{{\\dag}} = {latex(V_r)} {latex(S_r_inv)} {latex(U_r_T)}$\n")
    output_content.append(f"$A^{{\\dag}} = {latex(V_r * S_r_inv)} {latex(U_r_T)}$\n")
    output_content.append(f"$A^{{\\dag}} = {latex(A_dag)}$\n")

    # Guardar en archivo
    with open(archivo_salida, "w", encoding='utf-8') as file:
        file.writelines(output_content)

    return A_dag

# MP_minimos_cuadrados: Descripción no disponible.
def MP_minimos_cuadrados(A_dag, A, b, archivo_salida="solucion_minimos_cuadrados.txt"):   
    A = Matrix(A)
    b = Matrix(b)
    A_dag = Matrix(A_dag)
    
    output_content = []  # Lista para almacenar el contenido impreso
    
    display(Markdown(f"Buscamos resolver $Ax = b$:\n\n ${latex(A)} x = {latex(b)}$\n"))
    output_content.append(f"Buscamos resolver $Ax = b$:\n\n ${latex(A)} x = {vector_a_latex(b)}$\n")

    # Mostrar la inversa de Moore-Penrose
    display(Markdown(f"La inversa de Moore-Penrose es:\n$A^{{\\dag}} = {latex(A_dag)}$"))
    output_content.append(f"La inversa de Moore-Penrose es:\n$A^{{\\dag}} = {latex(A_dag)}$\n\n")

    # Calcular la solución de mínimos cuadrados
    x = A_dag * b
    display(Markdown("Luego, la solución de mínimos cuadrados de $Ax = b$ es:"))
    display(Markdown(f"$x = A^{{\\dag}} b$"))
    display(Markdown(f"$x = {latex(A_dag)} {latex(b)}$"))
    display(Markdown(f"$x = {latex(x)}$"))
    output_content.append("Luego, la solución de mínimos cuadrados de $Ax = b$ es:\n")
    output_content.append(f"$x = A^{{\\dag}} b$\n")
    output_content.append(f"$x = {latex(A_dag)} {vector_a_latex(b)}$\n")
    output_content.append(f"$x = {vector_a_latex(x)}$\n")

    # Guardar en archivo
    with open(archivo_salida, "w", encoding='utf-8') as file:
        file.writelines(output_content)

    return x


def inversa_moore_penroseRCOM(A, archivo_salida="inversa_moore_penrose.txt", rang = False, inv = False):
    A = Matrix(A)
    n, m = A.shape
    output_content = []  # Lista para almacenar el contenido impreso
    
    display(Markdown(f"Calcularemos la inversa de Moore-Penrose de \n\n $A = {latex(A)}$\n"))
    output_content.append(f"Calcularemos la inversa de Moore-Penrose de \n\n $A = {latex(A)}$\n")  
    
    # Determinar el rango de A
    rango = A.rank()
    if rang:
        display(Markdown(f"Apliacamos reduccion a $A^T$:"))
        output_content.append(f"Apliacamos reduccion a $A^T$:")
        output_content.append(f"${latex(A.T)}$:\n\n")
        reduccion_triangular(A.T)
        display(Markdown(f"De esta forma hay {rango} columnas LI de A, por lo que\n\n $\\varphi(A) = {rango}$."))
        output_content.append(f"De esta forma hay {rango} columnas LI de A, por lo que\n\n $\\varphi(A) = {rango}$.")
    else:    
        display(Markdown(f"El rango de la matriz $A$ es {rango}.\n"))
        output_content.append(f"El rango de la matriz $A$ es {rango}.\n\n")
        display(Markdown(f"$\\varphi(A) = {rango}$."))
        output_content.append(f"$\\varphi(A) = {rango}$.\n\n")

    # Calcular la inversa de Moore-Penrose según el rango
    if rango == m:
        # A^{\\dag} = (A^T A)^{-1} A^T
        
        A_dag = (A.T * A).inv() * A.T
        display(Markdown("De esta forma podemos usar:\n\n $ A^{\\dag} = (A^T A)^{-1} A^T $:\n"))
        output_content.append("Usando $A^{\\dag} = (A^T A)^{-1} A^T$:\n\n")
        display(Markdown(f"$A^{{\\dag}} = ({latex(A.T)} {latex(A)})^{{-1}} {latex(A.T)}$"))
        output_content.append(f"$A^{{\\dag}} = ({latex(A.T)} {latex(A)})^{{-1}} {latex(A.T)}$")
        display(Markdown(f"$A^{{\\dag}} = ({latex(A.T * A)})^{{-1}} {latex(A.T)}$"))
        output_content.append(f"$A^{{\\dag}} = ({latex(A.T * A)})^{{-1}} {latex(A.T)}$")
        if inv:
            display(Markdown("---"))
            display(Markdown("Calculamos la inversa de $A^T A$:\n"))
            output_content.append("Calculamos la inversa de $A^T A$:\n")
            inversa(A.T * A, archivo_salida="INVERSA_MP.txt")
            display(Markdown("---"))
        
        display(Markdown(f"$A^{{\\dag}} = {latex((A.T * A).inv())} {latex(A.T)}$"))
        output_content.append(f"$A^{{\\dag}} = {latex((A.T * A).inv())} {latex(A.T)}$")
        display(Markdown(f"$A^{{\\dag}} = {latex(A_dag)}$"))
        output_content.append(f"$A^{{\\dag}} = {latex(A_dag)}$")
        

    elif rango == n:
        # A^{\\dag} = A^T (A A^T)^{-1}
        A_dag = A.T * (A * A.T).inv()
        display(Markdown("De esta forma podemos usar:\n\n $A^{\\dag} = A^T (A A^T)^{-1}$:\n"))
        output_content.append("Usando $A^{\\dag} = A^T (A A^T)^{-1}$:\n\n")
        display(Markdown(f"$A^{{\\dag}} = {latex(A.T)} ({latex(A)} {latex(A.T)})^{{-1}}$"))
        output_content.append(f"$A^{{\\dag}} = {latex(A.T)} ({latex(A)} {latex(A.T)})^{{-1}}$")
        display(Markdown(f"$A^{{\\dag}} = {latex(A.T)} ({latex(A * A.T)})^{{-1}}$"))
        output_content.append(f"$A^{{\\dag}} = {latex(A.T)} ({latex(A * A.T)})^{{-1}}$")

        if inv:
            display(Markdown("---"))
            display(Markdown("Calculamos la inversa de $A A^T$:\n"))
            output_content.append("Calculamos la inversa de $A A^T$:\n")
            inversa(A * A.T, archivo_salida="INVERSA_MP.txt")
            display(Markdown("---"))

        display(Markdown(f"$A^{{\\dag}} = {latex((A * A.T).inv())} {latex(A.T)}$"))
        output_content.append(f"$A^{{\\dag}} = {latex((A * A.T).inv())} {latex(A.T)}$")
        display(Markdown(f"$A^{{\\dag}} = {latex(A_dag)}$"))
        output_content.append(f"$A^{{\\dag}} = {latex(A_dag)}$")


    else:
        raise ValueError("El rango de A no es igual a ni n ni m, no se puede calcular la inversa de Moore-Penrose con las fórmulas dadas.")

    # Mostrar y guardar la inversa de Moore-Penrose
    display(Markdown(f"La inversa de Moore-Penrose $A^{{\\dag}}$ es:\n\n$A^{{\\dag}} = {latex(A_dag)}$"))
    output_content.append(f"La inversa de Moore-Penrose $A^{{\\dag}}$ es:\n\n$A^{{\\dag}} = {latex(A_dag)}$\n")

    # Guardar en archivo
    with open(archivo_salida, "w", encoding='utf-8') as file:
        file.writelines(output_content)

    return A_dag

def vector_a_latex(v_i):
    """
    Convierte un vector columna de SymPy a una cadena de LaTeX representando el vector
    como una matriz columna.
    """
    # Comenzar la matriz en LaTeX
    latex_vector = "\\left[\\begin{matrix}"

    # Añadir cada elemento del vector a la matriz, en su propia línea
    for elemento in v_i:
        latex_vector += f" {latex(elemento)} \\\\"

    # Finalizar la matriz
    latex_vector += " \\end{matrix}\\right]"

    return latex_vector