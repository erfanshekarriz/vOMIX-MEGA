# Use an official Python runtime as a parent image
FROM python:3.11-slim

# Set the working directory in the container
WORKDIR /app

# Copy only the project metadata and source first for better caching
COPY pyproject.toml README.md /app/
COPY vomix /app/vomix

# Install build dependencies
RUN pip install --no-cache-dir --upgrade pip setuptools

# Install the package using PEP 517/518 standards
RUN pip install --no-cache-dir .

# Make port 8000 available to the world outside this container
EXPOSE 8000

# Define environment variable
ENV PYTHONUNBUFFERED=1

# Run the CLI entrypoint as defined in pyproject.toml
CMD ["vomix"]