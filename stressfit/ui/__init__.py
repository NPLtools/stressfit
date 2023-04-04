"""Package encapsulates user interface functionality."""
_nb = None
def run_gui():
    """Run Jupyter notebook GUI."""
    global _nb
    from .notebook import UI
    _nb = UI()
    _nb.display()

def _create_template(file='stressfit_template.py'):
    from pathlib import Path as _Path
    import shutil
    templ = file
    thisfn = _Path(__file__)
    p = thisfn.parent.parent
    src = p.joinpath('resources/scripts')
    src = src.joinpath(templ)
    tgt = src.cwd().joinpath(templ)
    shutil.copy(src, tgt)
    print('Template file copied in {}'.format(str(tgt)))

def create_notebook():
    """Export Jupyter notebook template file to the current directory."""
    _create_template(file='stressfit_template.ipynb')
    
def create_script():
    """Export python script template file to the current directory."""
    _create_template(file='stressfit_template.py')
    
