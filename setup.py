#!/usr/bin/env python
# file cubit/setup.py

"""Compilation script for the cubit package.

CUBIT - The CUBIT Utility BioInformatics Toolkit.  A collection of
bioinformatics scripts written by James Smagala for use by the CDC flu
branch.

Revision History:
    20081020 James Smagala: hacked together from many other setup scripts
    20090203 JS: Heavy revisions to replace make.bat
    20100701 JS: Remove cubit/cubit_extras distinction

"""

from distutils.core import setup
import py2exe
import sys
import os
import os.path
import subprocess
import shutil
import markdown
from datetime import datetime


def inno_script_setup(name, root, windows, install_path, icon,
                      version=datetime.now().strftime('%Y%m%d'),
                      iss_file_name=None):
    """Make an InnoScript compile script"""
    root = os.path.abspath(root)
    path = iss_file_name or name + '.iss'
    install_dir = os.path.join(install_path, name)

    # get a list of all the files to package
    all_files = []
    for walked_dir in os.walk(root):
        prefix = os.path.commonprefix((root, walked_dir[0], ))
        for file_name in walked_dir[2]:
            relative_path = walked_dir[0][len(prefix) + 1:]
            all_files.append(os.path.join(relative_path, file_name))

    # write out the script
    outfile = open(path, 'w')
    print >> outfile, '; WARNING: This script has been created by'
    print >> outfile, '; setup.py.  Changes to this script will be'
    print >> outfile, '; overwritten the next time setup.py is run!'
    print >> outfile, '[Setup]'
    print >> outfile, 'AppName=Install %s' % name
    print >> outfile, 'AppVerName=%s %s' % (name, version)
    print >> outfile, 'DefaultDirName=%s' % install_dir
    print >> outfile, 'DefaultGroupName=%s' % name
    print >> outfile, 'SetupIconFile=%s' % icon
    print >> outfile, ''

    print >> outfile, r"[Files]"
    for path in all_files:
        print >> outfile, r'Source: "%s"; DestDir: "{app}\%s"; Flags: '\
                          r'ignoreversion' % \
                          (os.path.join(name, path),
                           os.path.dirname(path))
    print >> outfile, ''

    print >> outfile, r"[Icons]"
    for window in windows:
        file_name = window['dest_base']
        path = file_name + '.exe'
        icon = os.path.normpath(window['icon_resources'][0][1])
        print >> outfile, r'Name: "{group}\%s"; Filename: '\
                          r'"{app}\%s"; IconFilename: "{app}\%s"' % \
                          (file_name, path, icon, )
        print >> outfile, r'Name: "{commondesktop}\%s"; Filename: '\
                          r'"{app}\%s"; WorkingDir: "{app}"; '\
                          r'IconFilename: "{app}\%s"' % \
                          (file_name, path, icon)
    print >> outfile, 'Name: "{group}\Uninstall %s"; Filename: '\
                      '"{uninstallexe}"; IconFilename: '\
                      '"{app}\icons\install.ico"'% name


def get_manifest():
    """Get a resouce manifest - gives apps a WinXP look if run on XP"""
    manifest_template = """
<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<assembly xmlns="urn:schemas-microsoft-com:asm.v1" manifestVersion="1.0">
  <assemblyIdentity
    version="5.0.0.0"
    processorArchitecture="x86"
    name="%(prog)s"
    type="win32"
  />
  <description>%(prog)s</description>
  <trustInfo xmlns="urn:schemas-microsoft-com:asm.v3">
    <security>
      <requestedPrivileges>
        <requestedExecutionLevel
            level="asInvoker"
            uiAccess="false">
        </requestedExecutionLevel>
      </requestedPrivileges>
    </security>
  </trustInfo>
  <dependency>
    <dependentAssembly>
      <assemblyIdentity
            type="win32"
            name="Microsoft.VC90.CRT"
            version="9.0.21022.8"
            processorArchitecture="x86"
            publicKeyToken="1fc8b3b9a1e18e3b">
      </assemblyIdentity>
    </dependentAssembly>
  </dependency>
  <dependency>
    <dependentAssembly>
        <assemblyIdentity
            type="win32"
            name="Microsoft.Windows.Common-Controls"
            version="6.0.0.0"
            processorArchitecture="X86"
            publicKeyToken="6595b64144ccf1df"
            language="*"
        />
    </dependentAssembly>
  </dependency>
</assembly>
"""
    rt_manifest = 24
    return rt_manifest, manifest_template


def get_windows():
    """Get defs for GUI windows"""
    rt_manifest, manifest_template = get_manifest()
    mutation_mapper = {
        'script': 'main/mutation_mapper.py',
        'icon_resources': [(1, 'icons/mutation_mapper.ico', ), ],
        'other_resources': [(rt_manifest, 1,
            manifest_template % {'prog': 'Mutation Mapper', }, ), ],
        'dest_base': 'Mutation Mapper',
    }
    aadiff = {
        'script': 'main/aadiff.py',
        'icon_resources': [(1, 'icons/aadiff.ico', ), ],
        'other_resources': [(rt_manifest, 1,
            manifest_template % {'prog': 'AA Diff', }, ), ],
        'dest_base': 'AA Diff',
    }
    count_gs = {
        'script': 'main/count_gs.py',
        'icon_resources': [(1, 'icons/count_gs.ico', ), ],
        'other_resources': [(rt_manifest, 1,
            manifest_template % {'prog': 'Count GS', }, ), ],
        'dest_base': 'Count GS',
    }
    fasta_tree_order = {
        'script': 'main/fasta_tree_order.py',
        'icon_resources': [(1, 'icons/fasta_tree_order.ico', ), ],
        'other_resources': [(rt_manifest, 1,
            manifest_template % {'prog': 'Fasta Tree Order', }, ), ],
        'dest_base': 'Fasta Tree Order',
    }
    swap_plates = {
        'script': 'main/swap_plates.py',
        'icon_resources': [(1, 'icons/swap_plates.ico', ), ],
        'other_resources': [(rt_manifest, 1,
            manifest_template % {'prog': 'Swap Plates', }, ), ],
        'dest_base': 'Swap Plates',
    }
    remap_fasta = {
        'script': 'main/remap_fasta.py',
        'icon_resources': [(1, 'icons/remap_fasta.ico', ), ],
        'other_resources': [(rt_manifest, 1,
            manifest_template % {'prog': 'Remap FASTA', }, ), ],
        'dest_base': 'Remap FASTA',
    }
    windows = [mutation_mapper, aadiff, count_gs, fasta_tree_order,
            remap_fasta, swap_plates]
    return windows


def get_console():
    """Get defs for console-only windows"""
    console = [os.path.join('bin', script) for script in os.listdir('bin')
				if not script.endswith('.svn')]
    return console


def get_options(project_name, includes):
    """Get the py2exe option dictionary"""
    dll_excludes = [
        'iconv.dll',
        'intl.dll',
        'libatk-1.0-0.dll',
        'libgdk_pixbuf-2.0-0.dll',
        'libgdk-win32-2.0-0.dll',
        'libglib-2.0-0.dll',
        'libgmodule-2.0-0.dll',
        'libgobject-2.0-0.dll',
        'libgthread-2.0-0.dll',
        'libgtk-win32-2.0-0.dll',
        'libpango-1.0-0.dll',
        'libpangowin32-1.0-0.dll',
        # comment this out if os.popen() is used
        'w9xpopen.exe',
        'msvcp90.dll',
        'msvcr90.dll',
    ]
    excludes = []
    packages = ['gui', 'main', ]

    # actually build the options dictionary
    options = {'py2exe': {
        'compressed': 1,
        'optimize': 2,
        'bundle_files': 1,
        #'ascii': 1, # this will exclude encodings, making
        #            # the file smaller, but it will
        #            # break if anything uses unicode
        'dist_dir': project_name,
        'dll_excludes': dll_excludes,
        'excludes': excludes,
        'includes': includes,
        'packages': packages,
    }, }

    return options


def get_data_files(project_name, includes):
    """Get any extra files that have not been included yet"""
    data_files = []

    # check if we need matplotlib
    if 'matplotlib' in includes:
        from matplotlib import get_py2exe_datafiles
        data_files.append(get_py2exe_datafiles())

    # include some MS dll files to reduce potential dependencies
    #dll_includes = ['misc/msvcp90.dll', 'misc/msvcr90.dll', ]
    dll_includes = ['dll/msvcp90.dll']
    data_files.append(('lib', dll_includes, ))

    # the zip file should be prebuilt in the build dir
    data_files.append(('source', ['build/%s_source.zip' % project_name, ], ))

    # any icons used by any of the programs go here
    icons = [os.path.join('icons', icon) for icon in os.listdir('icons') 
			if not icon.endswith('.svn')]
    data_files.append(('icons', icons, ))

    # docs are assumed to be html formatted
    docs = [os.path.join('docs', doc) for doc in os.listdir('docs')
            if doc.endswith('html') or doc.endswith('.png') or
	    doc.endswith('hlp')]
    data_files.append(('docs', docs, ))

    return data_files


def rm(path):
    """Remove a file, but don't error out if it doesn't exist"""
    try:
        os.remove(path)
    except OSError:
        pass


def rename(src, dst):
    """Remove a file, but don't error out if it doesn't exist"""
    try:
        os.rename(src, dst)
    except OSError:
        pass


def copy(src, dst):
    """Copy a file, but don't error out if it doesn't exist"""
    try:
        shutil.copy(src, dst)
    except OSError:
        pass
    except IOError:
        pass


def rmtree(path):
    """Remove a tree structure without error if it doesn't exist"""
    try:
        shutil.rmtree(path)
    except OSError:
        pass


def mkdir(path):
    """Create a directory, but don't error out if it already exists"""
    try:
        os.makedirs(path)
    except OSError:
        pass


def create_html_file(markdownfile, htmlfile):
    """Run markdown to convert to html, and add headers, etc"""
    outfile = open(htmlfile, 'w')
    infile = open(markdownfile, 'r')
    print >> outfile, '<?xml version="1.0"?>\n' + \
                      '<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 ' + \
                      'Strict//EN"\n' + \
                      ' "http://www.w3.org/TR/xhtml1/DTD/xhtml1-strict.dtd">'
    print >> outfile, '<html xmlns="http://www.w3.org/1999/xhtml">\n' + \
                      '<head>\n' + \
                      '<meta http-equiv="Content-Type" ' + \
                      'content="text/html; charset=utf-8" />\n' + \
                      '<title>%s</title>\n</head>\n<body>' % \
                      os.path.splitext(os.path.basename(markdownfile))[0]
    print >> outfile, markdown.markdown(infile.read())
    print >> outfile, '</body>\n</html>'

    infile.close()
    outfile.close()


def pre_compile_tasks():
    """Stuff to run before compiling"""

    def gen_html_docs():
        """Refresh the generated html documentation"""
        # regenerate the html documentation
        print 'Refresh the html documentation:',
        os.chdir('docs')
        for name in os.listdir(os.getcwd()):
            print '.',
            if name.endswith('.txt'):
                html_file = os.path.splitext(name)[0] + '.html'
                create_html_file(name, html_file)
        os.chdir('..')

    def build_cubit_source():
        """Copy over the source code, but don't copy revision history"""
        print '\nBuild the source .zip:',
        mkdir('build')
        top = os.getcwd()
        for root, dirs, files in os.walk(top):
            print '.',
            prefix = os.path.commonprefix((top, root, ))
            rel_root = root[len(prefix) + 1:]
            if root == top:
                for dir in ('build', 'cubit', '.svn'):
                    if dir in dirs:
                        del dirs[dirs.index(dir)]
            new_path = os.path.join('build', 'cubit', rel_root)
	    if new_path.find('.svn') > -1:
                continue
            mkdir(new_path)
            for name in files:
                if name.find('.svn') > -1 or name.startswith('.') or \
			name.endswith('.swp') or \
                        name.endswith('.pyc') or name.endswith('.html'):
                    continue
                shutil.copy(os.path.join(root, name),
                            os.path.join(new_path, name))

    def zip_source():
        """Zip the cubit source up so it can be included when compiling"""
        print '\nBuild the source.zip files:'
        os.chdir('build')
        subprocess.call(['zip', '-r', '-m', 'cubit_source.zip', 'cubit/', ])
        os.chdir('..')

    gen_html_docs()
    build_cubit_source()
    zip_source()


def post_compile_tasks(windows, project_name, install_path):
    """Package the compiled files for distribution"""
    # make an inno setup .iss script, compile it to .exe, and also zip
    # the files for manual installs
    root = project_name
    inno_script_setup(name=root, root=root,
                      windows=windows, install_path=install_path,
                      icon=r'icons\install.ico',
                      iss_file_name=project_name + '.iss')
    subprocess.call(['ISCC.exe',
                     r'/O%s' % os.getcwd(), r'/F%s' % project_name,
                     r'%s.iss' % project_name, ])
    subprocess.call(['zip', '-r', '-m', '%s.zip' % project_name, root, ])


def compile():
    """Run the setup/py2exe compilation"""
    # add the current dir to the python path
    sys.path.append(os.getcwd())

    # run the precompile tasks to get a source .zip file ready
    clean()
    pre_compile_tasks()

    # actually run the distutils setup
    print 'Compiling:'
    project_name = 'cubit'
    install_path = r'C:\BioEdit\apps'
    includes = ['wx', ]
    windows = get_windows()

    # as if we included py2exe on the command line - this will
    # allow us to run the script by double clicking
    sys.argv.append('py2exe')
    setup(console = get_console(),
            options = get_options(project_name, includes),
            windows = windows,
            zipfile = 'lib/%s.lib' % project_name,
            data_files = get_data_files(project_name, includes), )

    # run the post-compile tasks to move compiled files
    post_compile_tasks(windows, project_name, install_path)


def clean():
    """clean up previous compilations and build source bundles"""
    print '\nClean-up old files:',
    rmtree('cubit')
    print '.',
    rmtree('build')
    print '.',

    for root, dirs, files in os.walk(os.getcwd()):
        # clean up all other dirs of any generated files
        for name in files:
            if name.endswith('.html') or name.endswith('.iss') or \
                    name.endswith('.pyc') or name.endswith('.zip') or \
                    name.endswith('.exe') or name.endswith('~'):
                rm(os.path.join(root, name))
    print '.'


if __name__ == '__main__':

    def main():
        """Run the compliation if setup is double clicked or run from shell"""
        if 'clean' in sys.argv:
            clean()
        else:
            compile()

    main()
