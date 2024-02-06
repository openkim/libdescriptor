# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'Libdescriptor'
copyright = '2023, Amit Gupta'
author = 'Amit Gupta'
release = '0.2.0'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = ['sphinx.ext.autodoc', 'sphinx.ext.napoleon', 'sphinx.ext.viewcode', 'sphinx.ext.mathjax']

templates_path = ['_templates']
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']



# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'furo'
html_static_path = ['_static']
html_logo = '_static/libdescriptor_logo.png'
html_theme_options = {
    "sidebar_hide_name": True,
}
# -- Options for breathe -----------------------------------------------------
breath_projects = { "libdescriptor": "/home/amit/Projects/COLABFIT/colabfit-kim-model/colabfit-descriptor-library/tmp-doc/docs/xml" }
# TODO: change to relative path
