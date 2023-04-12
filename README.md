## Ejercicio final de transcriptómica (curso 2022-2023)
## Alumno: Marcos Rubio Fernández
## Abril 2023

A continuación se detalla resumidamente los pasos para este ejercico

### Apartado 1

Para llevar a cabo este apartado debe activarse en primer lugar el ambiente correspondiente de la carpeta envs

`conda env create -f Apartado1_env.yml`
`conda activate Apartado1_env`

El informe correspondiente se encuentra en la carpeta Memorias y contiene la explicación de cada parte del código(s) generados. Los diferentes scripts detallados en la memoria del ejercicio están localizados en el directorio scripts:

ejemplo:

`bash scripts/Apartado1_pipeline-basico.sh`

De esta manera se generarán de manera local todos los archivos necesarios. Algunos de los resultados se han subido al repositorio por si facilitan la evaluación del ejercicio. Para más información leer la memoria correspondiente.


### Apartado 2

Para la resolución de las dos preguntas correspondientes a este apartado se debe activar un único ambiente conda

`conda env create -f Apartado2_env.yml`
`conda activate Apartado2_env`

Para cada pregunta se han creado unos archivos .Rmd para la revisión del código:

- Pregunta 4: Apartado2_DEGtask.Rmd
- Pregunta 5: Apartado2_GSEA.Rmd

Además se proporciona en este directorio los archivos renderizados en .pdf y .html para la evaluación. También se han incluido en el repositorio los archivos correspondientes al resultado de GSEA.