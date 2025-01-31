from setuptools import setup, find_packages

setup(
    name="Constructor",  # Package name
    version="1.5",
    author="Martin Mihaylov",
    description="A description of the Constructor module",
    packages=find_packages(),  # This will include the 'Constructor' package
    install_requires=[],  # Add dependencies if needed
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires=">=3.6",
)