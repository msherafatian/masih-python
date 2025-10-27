# syntax=docker/dockerfile:1
ARG PYTHON_VERSION=3.11.13
FROM python:${PYTHON_VERSION}-slim AS base

ENV PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1

# System deps
RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    && rm -rf /var/lib/apt/lists/*

# Create non-root user
RUN useradd -m -u 1000 appuser

WORKDIR /app

# Copy and install requirements first (for better caching)
COPY requirements.txt .
RUN --mount=type=cache,target=/root/.cache/pip \
    python -m pip install --no-cache-dir -r requirements.txt

# Copy the entire application
COPY . .

# If you want to install masih as a package (optional, only if you have pyproject.toml)
# Uncomment the following line if you want the package installed in site-packages:
# RUN pip install --no-cache-dir -e .
# Create sessions directory with proper permissions
RUN mkdir -p /app/sessions && \
    chown -R appuser:appuser /app && \
    chmod 755 /app/sessions

# Switch to non-root user
USER appuser

EXPOSE 8050

# Run with gunicorn
CMD ["gunicorn", "app:server", "--bind=0.0.0.0:8050", "--workers=4", "--timeout=120"]