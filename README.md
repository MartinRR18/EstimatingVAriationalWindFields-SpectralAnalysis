# Estimating Variational Wind Fields via Spectral Analysis

## Descripción del Proyecto
Este repositorio contiene una investigación técnica y herramientas computacionales para estimar campos de viento ajustados mediante la resolución de ecuaciones elípticas. El enfoque principal es el uso de **Series de Fourier** para transformar problemas de variables continuas en problemas algebraicos manejables, permitiendo un análisis detallado de cómo las condiciones de contorno afectan la precisión en regiones de mesoescala.

El objetivo central es calcular el **porcentaje de pérdida de masa** en diferentes Problemas de Valor de Frontera (BVP, por sus siglas en inglés) para validar la robustez de los modelos de ajuste de viento.

## Contenido del Repositorio
El repositorio organiza sus archivos en cuadernos de Jupyter (`.ipynb`) y scripts de Python (`.py`) para facilitar tanto el análisis teórico como la ejecución numérica:

### Análisis de Problemas de Valor de Frontera (BVP)
Los archivos `BVP1_Ex3.ipynb`, `BVP2_Ex3.ipynb` y `BVP3_Ex3.ipynb` representan los experimentos principales:
*   **BVP1, BVP2, BVP3:** Cada uno explora diferentes configuraciones de contorno.
*   **Scripts de Soporte:** Se incluyen archivos `.py` específicos (`BVP1_Ex3.py`, `SpecialCaseNeumann.py`, `SpecialCaseMixtas.py`) que contienen las implementaciones algorítmicas de las soluciones para condiciones de tipo Neumann y condiciones mixtas.

### Herramientas de Cálculo y Visualización
*   **`CalculoCoeficientes.ipynb`**: Cuaderno dedicado al cálculo simbólico de los coeficientes de Fourier necesarios para la reconstrucción del campo de viento.
*   **`FourierMAtplotlib.ipynb`**: Notebook enfocado en la visualización gráfica de alta precisión de las series de Fourier resultantes.
*   **`SeriesFourierProy1.ipynb`**: Análisis fundamental sobre la convergencia y aproximación de las series.
*   **`Luf_CF_DMN.ipynb`**: Cuaderno especializado para el procesamiento de datos específicos (posiblemente relacionado con campos de viento de mesoescala).

## Tecnologías Utilizadas
El proyecto aprovecha un ecosistema robusto de bibliotecas científicas en Python:
*   **`SymPy`**: Para el cálculo simbólico y la derivación analítica de los coeficientes.
*   **`Matplotlib`**: Para la generación de gráficos detallados y el análisis de errores.
*   **`NumPy`**: Como motor principal para el cálculo numérico, generación de arreglos y operaciones matriciales eficientes.

## Instalación y Requisitos
Para ejecutar estos cuadernos, se recomienda tener instalado un entorno de Python con las siguientes dependencias:
*   Jupyter Notebook / JupyterLab
*   NumPy
*   SymPy
*   Matplotlib

## Uso Básico
1. Clonar el repositorio:
   ```bash
   git clone https://github.com/MartinRR18/EstimatingVAriationalWindFields-SpectralAnalysis.git

2. Instalar las dependencias (recomendado usar un entorno virtual):
   ```bash   
   pip install numpy sympy matplotlib jupyte
3.Abrir los cuadernos .ipynb con Jupyter Notebook o JupyterLab para explorar los análisis y ejecutar los cálculos.

### Metodología
El enfoque principal consiste en resolver una ecuación elíptica ajustando el campo de viento mediante series de Fourier. Se analizan tres diferentes Problemas de Valor de Frontera (BVP) con distintas condiciones de contorno para evaluar su impacto en la conservación de masa, calculando el porcentaje de pérdida de masa en cada caso.

Este análisis permite entender cómo las condiciones de frontera afectan la precisión y estabilidad de los modelos de viento en escalas meso, lo cual es crucial para aplicaciones meteorológicas y ambientales.

Si tienes alguna duda o quieres contribuir, no dudes en abrir un issue o enviar un pull request.

