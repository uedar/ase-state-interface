# ase-state-interface

## Installation

```
$ git clone git@github.com:uedar/ase-state-interface.git
$ cd ase-state-interface
$ pyenv local 3.7.5
$ poetry install
```

- [pyenv](https://uedar.github.io/install/pyenv.html)
- [poetry](https://uedar.github.io/install/poetry.html)

### Installation via Conda package manager
```
$ conda create -n ase_state_env python=3.7 pip -y
$ git clone https://github.com/uedar/ase-state-interface.git
$ conda activate ase_state_env
$ pip install ase-state-interface
```
---
## Developer Notes
### Setup developer environment using Conda

```bash
# Needs python3.6, ase for maximum compatability
# needs Proper linting tools (flake8 pydocstyle pylint autopep8) for compliance
# Create conda environment:
$ conda install -n dev -c conda-forge python=3.6 tk flask pylint flake8 ase=3.22 pytest=5.2 jupyterlab autopep8 pydocstyle
```

### Visual Studio Code (vscode) python linting configuration
- flake8, pylint and pydocstyle will be used together.
STEPS:
1. File > Preferences > Settings > Workspace
2. Search 'python'
3. Click 'Edit in settings.json'
4. Add the following code
```json
{
  "python.linting.flake8Enabled": true,
  "python.linting.pylintEnabled": true,
  "python.linting.enabled": true,
  "python.linting.pydocstyleEnabled": true
}
```
5. Save and restart


### Code-style Compliance checklist
- [ ] Run `pylint filename.py`
- [ ] Run `flake8 filename.py`