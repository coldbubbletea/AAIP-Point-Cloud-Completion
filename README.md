## Work for submission of ICDE24







### Usage

#### 1) Envrionment & prerequisites

- Pytorch 1.12.1
- CUDA 11.6
- Python 3.10
- [Visdom](https://github.com/facebookresearch/visdom)
- [Open3D](http://www.open3d.org/docs/release/index.html#python-api-index)

#### 2) Compile

Compile our EMD modules:  

    cd emd
    python3 setup.py install

---
header-includes:
  - \usepackage[ruled,vlined,linesnumbered]{algorithm2e}
---
# Algorithm 1
Just a sample algorithmn
\begin{algorithm}[H]
\DontPrintSemicolon
\SetAlgoLined
\KwResult{Write here the result}
\SetKwInOut{Input}{Input}\SetKwInOut{Output}{Output}
\Input{Write here the input}
\Output{Write here the output}
\BlankLine
\While{While condition}{
    instructions\;
    \eIf{condition}{
        instructions1\;
        instructions2\;
    }{
        instructions3\;
    }
}
\caption{While loop with If/Else condition}
\end{algorithm} 
