e=1
alpha=30
all: trayectoria_xy_$(e)_$(alpha).png trayectoria3D_$(e)_$(alpha).png particle_in_field.x graficas.py trayectoria_$(e)_$(alpha).dat

particle_in_field.x:particle_in_field.c graficas.py
	echo "Generando el archivo ejecutable"
	touch particle_in_field.x	
	cc particle_in_field.c -o particle_in_field.x
trayectoria_$(e)_$(alpha).dat:particle_in_field.x
	echo "Creando archivo de datos"
	./particle_in_field.x $(e) $(alpha)
trayectoria_xy_$(e)_$(alpha).png:trayectoria_$(e)_$(alpha).dat
	echo "Creando grafica 2D"
	python graficas.py trayectoria_$(e)_$(alpha).dat
trayectoria3D_$(e)_$(alpha).png:trayectoria_$(e)_$(alpha).dat
	echo "Creando grafica 3D"
	python graficas.py trayectoria_$(e)_$(alpha).dat
