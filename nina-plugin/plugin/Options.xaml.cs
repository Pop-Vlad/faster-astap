using System.ComponentModel.Composition;
using System.Windows;

namespace FasterAstap {

    /// <summary>
    /// Exporting the dictionary is what makes N.I.N.A. merge it and find the
    /// options DataTemplate inside.
    /// </summary>
    [Export(typeof(ResourceDictionary))]
    public partial class Options : ResourceDictionary {

        public Options() {
            InitializeComponent();
        }
    }
}
