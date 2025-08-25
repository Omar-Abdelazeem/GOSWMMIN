# GOSWMMIN Installation Troubleshooting Guide

This guide helps resolve common installation issues with the GOSWMMIN package.

## Common Issues and Solutions

### 1. Conda Environment Creation Fails

#### Problem
```bash
conda env create -f GOSWMMIN.yaml
# Error: PackagesNotFoundError or ResolvePackageNotFound
```

#### Cause
The original `GOSWMMIN.yaml` file has several issues:
- Platform-specific build strings (Windows-only packages on macOS/Linux)
- Overly restrictive version constraints
- Outdated package versions
- Missing channels or conflicts

#### Solution
Use the updated environment files:

```bash
# Use the new cross-platform environment
conda env create -f environment.yaml

# Or use the minimal version
conda env create -f environment-minimal.yaml
```

### 2. Package Version Conflicts

#### Problem
```bash
# Conflict resolution messages or dependency errors
```

#### Cause
Different packages requiring incompatible versions of dependencies.

#### Solution A: Use mamba (faster solver)
```bash
conda install mamba -c conda-forge
mamba env create -f environment.yaml
```

#### Solution B: Manual step-by-step installation
```bash
conda create -n goswmmin python=3.9
conda activate goswmmin
conda install -c conda-forge numpy pandas matplotlib tqdm jupyter
pip install wntr pyswmm
```

### 3. WNTR Installation Issues

#### Problem
```bash
# WNTR fails to install or has dependency conflicts
```

#### Cause
WNTR has specific requirements for NetworkX and other graph libraries.

#### Solution
```bash
# Install WNTR dependencies first
conda install -c conda-forge networkx scipy
pip install wntr

# Or use specific versions if needed
pip install "wntr>=1.0.0,<2.0"
```

### 4. PySWMM Installation Issues

#### Problem
```bash
# PySWMM installation fails or SWMM not found
```

#### Cause
- Missing system-level SWMM installation
- Compiler issues on some systems

#### Solution A: Use conda-forge version
```bash
conda install -c conda-forge pyswmm
```

#### Solution B: Install SWMM system-wide first
- **Windows**: Download SWMM from EPA website
- **macOS**: `brew install swmm`
- **Linux**: Install from package manager or compile from source

### 5. Platform-Specific Issues

#### Windows
```bash
# If you encounter M2W64 errors, create fresh environment
conda create -n goswmmin python=3.9
conda activate goswmmin
conda install -c conda-forge numpy pandas matplotlib pyswmm
pip install wntr tqdm
```

#### macOS
```bash
# For Apple Silicon Macs, ensure you're using the right architecture
conda create -n goswmmin python=3.9
conda activate goswmmin
conda install -c conda-forge numpy pandas matplotlib
pip install wntr pyswmm tqdm
```

#### Linux
```bash
# Install build essentials if needed
sudo apt-get update
sudo apt-get install build-essential python3-dev

# Then proceed with conda environment
conda env create -f environment.yaml
```

### 6. Import Errors After Installation

#### Problem
```python
import goswmmin
# ImportError or ModuleNotFoundError
```

#### Cause
- Package not installed in development mode
- Wrong Python environment activated

#### Solution
```bash
# Ensure you're in the right environment
conda activate goswmmin

# Install in development mode
cd /path/to/GOSWMMIN/pkg
pip install -e .

# Verify installation
python -c "import goswmmin; print(goswmmin.__version__)"
```

### 7. Jupyter Kernel Issues

#### Problem
```bash
# Jupyter doesn't see the conda environment
```

#### Solution
```bash
conda activate goswmmin
conda install ipykernel
python -m ipykernel install --user --name goswmmin --display-name "GOSWMMIN"
```

## Step-by-Step Fresh Installation

If all else fails, here's a complete fresh installation process:

### Step 1: Clean Installation
```bash
# Remove any existing environments
conda env remove -n goswmmin
conda env remove -n GOSWMMIN

# Clear conda cache
conda clean -a
```

### Step 2: Create New Environment
```bash
# Create minimal environment
conda create -n goswmmin python=3.9 -y
conda activate goswmmin
```

### Step 3: Install Core Dependencies
```bash
# Install via conda-forge (more up-to-date)
conda install -c conda-forge numpy pandas matplotlib tqdm jupyter -y
```

### Step 4: Install Specialized Packages
```bash
# Install WNTR and PySWMM via pip
pip install wntr>=1.0.0
pip install pyswmm>=2.0.0
```

### Step 5: Install GOSWMMIN
```bash
cd /path/to/GOSWMMIN/pkg
pip install -e .
```

### Step 6: Verify Installation
```bash
python -c "
import goswmmin
import wntr
import pyswmm
print('✓ All packages imported successfully')
print(f'GOSWMMIN version: {goswmmin.__version__}')
"
```

## Alternative: Docker Installation

For a completely isolated installation:

```dockerfile
FROM continuumio/miniconda3

WORKDIR /app

# Copy environment file
COPY environment-minimal.yaml .

# Create environment
RUN conda env create -f environment-minimal.yaml

# Activate environment
SHELL ["conda", "run", "-n", "goswmmin-minimal", "/bin/bash", "-c"]

# Copy and install package
COPY . .
RUN pip install -e .

# Set default command
CMD ["conda", "run", "-n", "goswmmin-minimal", "python", "-c", "import goswmmin; print('GOSWMMIN ready!')"]
```

## Getting Help

If you continue to experience issues:

1. **Check the GitHub Issues**: [GOSWMMIN Issues](https://github.com/omaraliamer/GOSWMMIN/issues)
2. **Create a new issue** with:
   - Your operating system and version
   - Python version (`python --version`)
   - Conda version (`conda --version`)
   - Complete error message
   - Steps you've already tried

3. **Contact**: omaraliamer98@gmail.com

## Environment File Comparison

| File | Purpose | Size | Dependencies |
|------|---------|------|--------------|
| `GOSWMMIN.yaml` | Original (problematic) | Large | Windows-specific, overly constrained |
| `environment.yaml` | Recommended full setup | Medium | Cross-platform, development tools |
| `environment-minimal.yaml` | Minimal installation | Small | Just core dependencies |

**Recommendation**: Use `environment.yaml` for development work, `environment-minimal.yaml` for production deployments.
