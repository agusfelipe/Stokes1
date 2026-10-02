# Stokes1

Solver de elementos finitos en MATLAB para el problema de la cavidad con tapa móvil (*lid-driven cavity*) en el dominio [0,1]×[0,1]. La pared superior se mueve con velocidad horizontal 1; las otras tres paredes están en reposo (condición de no deslizamiento).

## Programas principales

- `mainStokes.m`: flujo de Stokes (lineal), −μΔu + ∇p = f, ∇·u = 0, con μ = 1.
- `mainNavierStokes_Picard.m`: Navier–Stokes estacionario (Re = 100) con linealización de Picard. Itera hasta que el incremento de velocidad y el residuo sean menores que 0.5e-8, con un máximo de 100 iteraciones. Pide opcionalmente un `.mat` con la velocidad inicial; si se deja vacío, arranca desde cero con la tapa en 1.

> **Pendiente:** `mainNavierStokes_Picard.m` llama a `ConvectionMatrix`, que todavía no está en el repositorio, por lo que no se puede ejecutar.

## Uso

Desde MATLAB, en la carpeta del proyecto:

```matlab
mainStokes
```

Pregunta el número de elementos por dirección (por defecto 10) y muestra las mallas, el campo de velocidades, las líneas de corriente y la presión.

## Elementos disponibles

Se configuran al inicio de cada `main*.m` mediante `elemV`, `degreeV` y `degreeP`:

| `elemV` | Elemento | Velocidad | Presión |
|---|---|---|---|
| 0 | Cuadriláteros | Q1, Q2 | Q1 (o constante con `degreeP=0`) |
| 1 | Triángulos | P1, P2 | P1 (o constante con `degreeP=0`) |
| 11 | Triángulos con burbuja (MINI) | P1 + burbuja | P1 |

Combinaciones estables recomendadas: Q2-Q1 (`elemV=0, degreeV=2, degreeP=1`) o MINI (`elemV=11, degreeV=1, degreeP=1`). La configuración por defecto de `mainStokes.m` (P1-P1) no cumple la condición inf-sup y produce oscilaciones en la presión.

## Estructura

| Archivo | Función |
|---|---|
| `CreateMeshes.m` | Genera las mallas de velocidad y presión |
| `CreateUniformMesh.m` | Malla estructurada uniforme |
| `CreateAdaptedMesh.m` | Malla refinada hacia las paredes (desactivada por defecto) |
| `StokesSystem.m` | Ensambla K (laplaciano), G (divergencia) y f (fuerzas másicas) |
| `SourceTerm.m` | Fuerza volumétrica (nula por defecto) |
| `BC_red.m` | Condiciones Dirichlet por reducción del sistema |
| `BC.m` | Alternativa con multiplicadores de Lagrange (no usada) |
| `Func_ReferenceElement/` | Elemento de referencia, cuadratura y funciones de forma |
| `PlotMesh.m`, `PlotResults.m`, `PlotStreamlines.m` | Visualización |
| `cinput.m` | Entrada por teclado con valor por defecto |

## Notas

- El flujo es confinado: la presión queda definida salvo una constante, que se fija en el nodo 1 (esquina inferior izquierda) con p = 0.
- El sistema de punto silla `[K Gᵀ; G 0]` se resuelve con un método directo (`A\b`).
- `PlotStreamlines.m` usa `tri2grid`, del PDE Toolbox de MATLAB.
