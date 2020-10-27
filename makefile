all: fluid.so
	@echo "Running"
	@python3 run.py

fluid.so: fluid.c
	@echo "Compiling"
	@cc -fPIC -shared -o fluid.so fluid.c