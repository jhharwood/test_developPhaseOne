FROM jupyter/datascience-notebook

RUN pip install --no-cache \
	scipy \
	numpy \
	pandas \
	pyproj \
	easygui \
	folium \
	affine \
	PySimpleGUI \
	utm 

COPY . /work

CMD ["jupyter", "notebook", "--port=8888", "--no-browser", "--ip=0.0.0.0", "--allow-root"]


