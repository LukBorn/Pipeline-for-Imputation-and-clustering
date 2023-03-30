import pip
package_names=['numpy', 'pandas', 'tk'] #packages to install
pip.main(['install'] + package_names + ['--upgrade'])