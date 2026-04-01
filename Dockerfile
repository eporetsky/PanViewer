FROM continuumio/miniconda3

WORKDIR /app

# Copy environment spec and create the conda env used by the app
COPY environment.yml .
RUN conda env create -f environment.yml && conda clean -afy

# Use the conda env from environment.yml (name: panviewer)
ENV PATH=/opt/conda/envs/panviewer/bin:$PATH
ENV CONDA_DEFAULT_ENV=panviewer

# Install a production WSGI server for Flask
RUN pip install --no-cache-dir gunicorn

# Application + data layout: panbarley.db should live at input/barley/panbarley.db
# (run build_index.py before docker build, or bake COPY of that path here).
COPY . /app

# The Flask app will be served by Gunicorn on port 80
EXPOSE 80
# app:app refers to "app" (module) and "app" (Flask instance) in app.py
CMD ["gunicorn", "-w", "4", "-b", "0.0.0.0:80", "app:app"]