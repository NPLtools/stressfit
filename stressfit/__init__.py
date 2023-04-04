from .dataio import __package_info__
__version__ = __package_info__['version']
__license__ = __package_info__['license']

def gui():
    """Run Jupyter notebook GUI."""
    from .ui import run_gui
    run_gui()
    
def create_notebook():
    """Export Jupyter notebook template file to the current directory."""
    from .ui import create_notebook
    create_notebook()
    
def create_script():
    """Export Jupyter notebook template file to the current directory."""
    from .ui import create_script
    create_script()
    
