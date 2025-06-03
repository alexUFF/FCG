from modules.Vetor import Vetor
import matplotlib.pyplot as plt

class Circulo:
    def __init__(self, centro, raio):
        if not isinstance(centro, Vetor):
            raise TypeError("O centro deve ser um objeto da classe Vetor.")
        if centro.dim != 2:
            raise ValueError("O círculo só pode ser definido em ℝ².")
        if not isinstance(raio, (int, float)) or raio <= 0:
            raise ValueError("O raio deve ser um número positivo.")

        self.centro = centro
        self.raio = raio

    def __str__(self):
        return f"Circulo de centro ({self.centro.x}, {self.centro.y}) e raio {self.raio}"

    def plot(self, fig=None, cor='purple', mostrar_eixos=True, limite=None):
        if fig is None:
            fig = plt.figure()
            ax = fig.add_subplot(1, 1, 1)
        else:
            ax = fig.gca()

        circulo = plt.Circle((self.centro.x, self.centro.y), self.raio,
                             edgecolor=cor, facecolor='none', linewidth=1.5, label=str(self))
        ax.add_patch(circulo)
        ax.plot(self.centro.x, self.centro.y, 'o', color=cor)  # marca o centro

        r = self.raio if limite is None else limite
        ax.set_aspect('equal')
        ax.set_xlim(self.centro.x - r*1.5, self.centro.x + r*1.5)
        ax.set_ylim(self.centro.y - r*1.5, self.centro.y + r*1.5)

        if mostrar_eixos:
            ax.axhline(0, color='k', linewidth=0.5)
            ax.axvline(0, color='k', linewidth=0.5)

        ax.legend()
        return fig
