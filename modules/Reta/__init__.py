from modules.Vetor import Vetor
import math

class Reta:
    def __init__(self, ponto, diretor=None, normal=None):
        if not isinstance(ponto, Vetor):
            raise TypeError("O ponto deve ser um objeto da classe Vetor.")
        if ponto.dim != 2 and diretor is None:
            raise ValueError("Retas 3d precisam de um vetor diretor")
        
        if (diretor is None and normal is None) or (diretor is not None and normal is not None):
            raise ValueError("Deve ser fornecido exatamente um vetor: diretor ou normal.")

        if diretor is not None:
            if not isinstance(diretor, Vetor) or diretor.dim != ponto.dim:
                raise ValueError("O vetor diretor tem que ter a mesma dimensão do ponto.")
            self.diretor = diretor
            if ponto.dim == 2:
                self.normal = Vetor([-diretor.y, diretor.x])  # Ortogonal ao vetor diretor
        else:
            if not isinstance(normal, Vetor) or normal.dim != 2:
                raise ValueError("O vetor normal deve ser um Vetor 2D.")
            # O vetor diretor é ortogonal ao vetor normal
            self.diretor = Vetor([-normal.y, normal.x])
            self.normal = normal

        self.ponto = ponto
        self.dim = ponto.dim

    def __str__(self, tipo='cartesiana'):
        x0, y0, z0  = self.ponto.x, self.ponto.y, self.ponto.z if self.ponto.dim == 3 else (0, 0)
        dx, dy, dz = self.diretor.x, self.diretor.y, self.diretor.z if self.diretor.dim == 3 else (0, 0)

        if tipo == 'cartesiana':
            if self.dim !=2 :
                raise ValueError("A reta cartesiana só pode ser definida em 2D.")
            
            return f"Reta: {self.normal.x}x + {self.normal.y}y = {self.normal.produto_interno(self.ponto)}"
        if tipo == "parametrica":  # tipo == 'normal'
            if self.dim == 2:
                return f"Reta: ({self.ponto.x} + t * {self.diretor.x}, {self.ponto.y} + t * {self.diretor.y})" 
            elif self.dim == 3:
                return f"Reta: ({self.ponto.x} + t * {self.diretor.x}, {self.ponto.y} + t * {self.diretor.y}, {self.ponto.z} + t * {self.diretor.z})"


    def plot(self, fig=None, cor='b', limite=5, mostrar_eixos=True):
        import matplotlib.pyplot as plt

        if fig is None:
            fig = plt.figure()
            ax = fig.add_subplot(1, 1, 1)
        else:
            ax = fig.gca()

        x0, y0 = self.ponto.x, self.ponto.y
        dx, dy = self.diretor.x, self.diretor.y

        t_vals = [-limite, limite]
        x_vals = [x0 + t * dx for t in t_vals]
        y_vals = [y0 + t * dy for t in t_vals]

        ax.plot(x_vals, y_vals, color=cor, label=str(self))
        ax.plot(x0, y0, 'o', color=cor)

        if mostrar_eixos:
            ax.axhline(0, color='k', linewidth=0.5)
            ax.axvline(0, color='k', linewidth=0.5)
            ax.set_aspect('equal')
            ax.set_xlim(x0 - limite, x0 + limite)
            ax.set_ylim(y0 - limite, y0 + limite)

        return fig

    def interseccao(self, outro):
        from modules.Circulo import Circulo
        from modules.Vetor import Vetor

        if self.dim != 2:
            raise ValueError("A interseção só está implementada para objetos no plano (dimensão 2).")

        # Interseção com outra reta
        if isinstance(outro, Reta):
            a1, b1 = self.normal.x, self.normal.y
            c1 = self.normal.produto_interno(self.ponto)

            a2, b2 = outro.normal.x, outro.normal.y
            c2 = outro.normal.produto_interno(outro.ponto)

            det = a1 * b2 - a2 * b1
            if det == 0:
                raise ValueError("As retas são paralelas ou coincidentes (sem interseção única).")

            x = (c1 * b2 - c2 * b1) / det
            y = (a1 * c2 - a2 * c1) / det
            return Vetor([x, y])

        # Interseção com círculo
        elif isinstance(outro, Circulo):
            # Parametrização da reta: X(t) = P + tV
            px, py = self.ponto.x, self.ponto.y
            dx, dy = self.diretor.x, self.diretor.y

            cx, cy = outro.centro.x, outro.centro.y
            r = outro.raio

            # Substitui na equação do círculo: ||X(t) - C||^2 = r^2
            a = dx**2 + dy**2
            b = 2 * (dx * (px - cx) + dy * (py - cy))
            c = (px - cx)**2 + (py - cy)**2 - r**2

            delta = b**2 - 4 * a * c

            if delta < 0:
                return []  # Sem interseção real
            elif abs(delta) < 10**(-9):
                t = -b / (2 * a)
                x = px + t * dx
                y = py + t * dy
                return [Vetor([x, y])]
            else:
                sqrt_delta = math.sqrt(delta)
                t1 = (-b + sqrt_delta) / (2 * a)
                t2 = (-b - sqrt_delta) / (2 * a)
                p1 = Vetor([px + t1 * dx, py + t1 * dy])
                p2 = Vetor([px + t2 * dx, py + t2 * dy])
                return [p1, p2]

        else:
            raise TypeError("interseccao espera outro objeto da classe Reta ou Circulo.")
